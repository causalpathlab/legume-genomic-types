//! Scoring a tool's cells, binning them, and the model that bins new cells.

use super::{compose, plain_batch, ColumnBlocks, CovariateArgs, Selection};
use crate::gene_sets::{CovariateSets, GeneSet, SetKind};
use crate::name_matching::GeneIndex;
use legume_numeric::matrix::covariate::{
    fit_bins, functionals, BinOpts, CovariateBins, Functional, GeneMeans, ScoreWeights, SetRows,
};
use log::info;
use nalgebra::DMatrix;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};

/// A set is scored when this many of its genes match, or half of a smaller
/// set (one gene of a one-gene set, XIST).
pub const MIN_MATCHED: usize = 10;
/// About the most cells the genes' mean expression is taken over.
const MEAN_CELLS: usize = 50_000;
/// The seed of the controls' draw and the mixture's starts.
const SEED: u64 = 42;

/// What a fit's covariates were and how its cells were binned, recorded
/// beside its run.
#[derive(Debug, Serialize, Deserialize)]
pub struct CovariateModel {
    /// Every set scored, by its label, and how (weights by gene name).
    pub functionals: Vec<Functional>,
    /// How the scores were binned.
    pub bins: CovariateBins,
    /// Cells per composite label.
    pub levels: BTreeMap<Box<str>, usize>,
    /// The names of the rows scored (a cell's totals come from them); new
    /// cells are scored on their rows of these names. Empty: every row.
    #[serde(default)]
    pub rows: Vec<Box<str>>,
}

/// The covariates as modelled: each cell's `batch|cvK` label in cell order,
/// every cell's scores (cells × functionals), and the model.
pub struct Modelled {
    pub labels: Vec<Box<str>>,
    pub scores: DMatrix<f32>,
    pub model: CovariateModel,
}

/// Each cell of `data`'s `batch|cvK` label and the fitted model; `None`
/// (the reason logged) when nothing is modelled: `none`, no covariate's
/// genes match, none with groups of its own, or one bin.
///
/// - `batch`: each cell's batch, in column order (a label an earlier fit
///   made `batch|cvK` counts as its plain batch).
/// - `proj_kn`: the cells' random projection (dims × cells) under those
///   labels, for the groups of like cells.
/// - `scored`: the rows scored and matched to gene sets (a cell's totals
///   and the control pool come from them only); all rows when `None`.
/// - `block`: columns read at a time.
///
/// # Errors
/// A malformed `--covariate`, an unknown covariate, a gene table or block
/// that cannot be read, or inputs whose sizes disagree.
pub fn model_covariates(
    data: &dyn ColumnBlocks,
    batch: &[Box<str>],
    proj_kn: &DMatrix<f32>,
    args: &CovariateArgs,
    scored: Option<&[usize]>,
    block: usize,
) -> anyhow::Result<Option<Modelled>> {
    let selection = args.selection()?;
    if matches!(selection, Selection::None) {
        return Ok(None);
    }
    anyhow::ensure!(
        batch.len() == data.n_cols(),
        "--covariate: {} batch labels for {} cells",
        batch.len(),
        data.n_cols()
    );
    let named = matches!(selection, Selection::Named(_));
    let batch = match plain_batches(batch) {
        Ok(b) => b,
        // Labels that already have parts of their own: `all` leaves them be.
        Err(e) if !named => {
            info!("--covariate all: not modelled: {e}");
            return Ok(None);
        }
        Err(e) => return Err(e),
    };
    let table = gene_sets(args)?;
    let chosen: Vec<&str> = match &selection {
        Selection::Named(names) => {
            let known = table.covariates();
            for n in names {
                anyhow::ensure!(
                    known.contains(&n.as_ref()),
                    "--covariate: no covariate {n} (known: {})",
                    known.join(", ")
                );
            }
            names.iter().map(AsRef::as_ref).collect()
        }
        _ => table.covariates(),
    };

    // Genes matched among the scored rows only; a covariate is scored when
    // every one of its sets matches.
    let names = data.row_names()?;
    let lookup = RowLookup::new(&names, scored);
    let (rows, index) = (&lookup.rows, &lookup.index);
    let mut sets: Vec<(&GeneSet, Vec<usize>)> = Vec::new();
    for c in &chosen {
        let found: Vec<_> = table.of(c).map(|s| (s, matched(s, index))).collect();
        if found.iter().all(|(_, r)| r.is_some()) {
            for (s, r) in found {
                let r: Vec<usize> = r.expect("matched").iter().map(|&i| rows[i]).collect();
                sets.push((s, r));
            }
        } else {
            info!("--covariate: {c} not scored, too few of its genes among the features");
        }
    }
    if sets.is_empty() {
        info!("--covariate: not modelled, no covariate's genes match the features");
        return Ok(None);
    }

    let means = if sets.iter().any(|(s, _)| s.kind == SetKind::Module) {
        Some(gene_means(data, names.len(), scored, block)?)
    } else {
        None
    };
    let labels_of: Vec<String> = sets.iter().map(|(s, _)| s.label()).collect();
    let set_rows: Vec<SetRows<'_>> = sets
        .iter()
        .zip(&labels_of)
        .map(|((s, r), label)| SetRows {
            label,
            kind: s.kind,
            rows: r,
        })
        .collect();
    let functionals = functionals(&set_rows, &names, means.as_deref(), scored, SEED);
    let scores = score(data, &functionals, &lookup, scored, block)?;
    let covariate_of: Vec<Box<str>> = sets.iter().map(|(s, _)| s.covariate.clone()).collect();
    let opts = BinOpts {
        judge: !named,
        fixed_k: args.covariate_bins,
        k_max: args.covariate_bins_max,
        criterion: args.covariate_criterion,
        seed: SEED,
        ..BinOpts::default()
    };
    let (bin, bins) = match fit_bins(&scores, &covariate_of, &batch, proj_kn, &opts)? {
        Ok(fit) => fit,
        Err(why) => {
            info!("--covariate: not modelled, {why}");
            return Ok(None);
        }
    };
    if !bins.set_aside.is_empty() {
        info!(
            "--covariate all: set aside {} (no groups of their own)",
            bins.set_aside.join(", ")
        );
    }

    let labels: Vec<Box<str>> = batch
        .iter()
        .zip(&bin)
        .map(|(b, &c)| compose(b, c))
        .collect::<anyhow::Result<_>>()?;
    let levels = labels.iter().fold(BTreeMap::new(), |mut m, l| {
        *m.entry(l.clone()).or_default() += 1;
        m
    });
    info!(
        "--covariate {}: {} bins over {} (own bins by covariate: {:?}; {:?} by k: {:?}); {} cells \
         kept plain batch labels ({} group(s) at least {:.0}% one bin)",
        chosen.join(","),
        bins.mixture.k(),
        bins.binned
            .iter()
            .map(|&j| functionals[j].label.as_ref())
            .collect::<Vec<_>>()
            .join(", "),
        bins.own_k,
        opts.criterion,
        bins.bic,
        bins.plain_cells,
        bins.one_bin_groups,
        100.0 * opts.one_bin
    );
    Ok(Some(Modelled {
        labels,
        scores,
        model: CovariateModel {
            functionals,
            bins,
            levels,
            rows: rows.iter().map(|&r| names[r].clone()).collect(),
        },
    }))
}

impl CovariateModel {
    /// Each cell of `data`'s `batch|cvK` label, scored and binned as the
    /// fit's cells were (a cell of a group the fit kept plain is binned: new
    /// cells form no groups). A batch the fit saw is centred on its recorded
    /// median and follows its merges; a batch it did not see is centred on
    /// its own median, so a small unseen batch gets a noisy centre, and takes
    /// its nearest bins as they are. `scored` and `block` as in
    /// [`model_covariates`]; without `scored`, the rows named as the fit's
    /// scored rows ([`Self::rows`]).
    ///
    /// # Errors
    /// A block that cannot be read, or a label that cannot be composed.
    pub fn labels(
        &self,
        data: &dyn ColumnBlocks,
        batch: &[Box<str>],
        scored: Option<&[usize]>,
        block: usize,
    ) -> anyhow::Result<Vec<Box<str>>> {
        let plain = plain_batches(batch)?;
        let names = data.row_names()?;
        let recorded: Option<Vec<usize>> = (scored.is_none() && !self.rows.is_empty()).then(|| {
            let fitted: std::collections::HashSet<&str> =
                self.rows.iter().map(AsRef::as_ref).collect();
            (0..names.len())
                .filter(|&r| fitted.contains(names[r].as_ref()))
                .collect()
        });
        let scored = scored.or(recorded.as_deref());
        let scores = score(
            data,
            &self.functionals,
            &RowLookup::new(&names, scored),
            scored,
            block,
        )?;
        self.bins
            .assign(&scores, &plain)
            .into_iter()
            .zip(&plain)
            .map(|(bin, b)| compose(b, bin))
            .collect()
    }

    /// The model a run recorded at `path`, if there is one.
    ///
    /// # Errors
    /// A recorded model that cannot be read.
    pub fn read(path: &std::path::Path) -> anyhow::Result<Option<Self>> {
        if !path.is_file() {
            return Ok(None);
        }
        let text = std::fs::read_to_string(path)?;
        Ok(Some(serde_json::from_str(&text).map_err(|e| {
            anyhow::anyhow!("{}: not a covariate model this reads: {e}", path.display())
        })?))
    }

    /// Record the model at `path`, as JSON.
    ///
    /// # Errors
    /// A file that cannot be written.
    pub fn write(&self, path: &std::path::Path) -> anyhow::Result<()> {
        std::fs::write(path, serde_json::to_string_pretty(self)?)?;
        Ok(())
    }
}

/// The built-in gene sets with `--covariate-genes`'s.
fn gene_sets(args: &CovariateArgs) -> anyhow::Result<CovariateSets> {
    let builtin = CovariateSets::builtin();
    Ok(match args.covariate_genes.as_deref() {
        Some(p) => {
            builtin.with(CovariateSets::read(std::path::Path::new(p)).map_err(anyhow::Error::msg)?)
        }
        None => builtin,
    })
}

/// The rows of `set` among `index`'s, when enough of its genes match.
fn matched(set: &GeneSet, index: &GeneIndex) -> Option<Vec<usize>> {
    let rows = set.matched_rows(index);
    let need = MIN_MATCHED.min(set.genes.len().div_ceil(2)).max(1);
    (rows.len() >= need).then_some(rows)
}

/// Each label without the covariate axis an earlier fit gave it, each
/// distinct label worked out once.
fn plain_batches(batch: &[Box<str>]) -> anyhow::Result<Vec<Box<str>>> {
    let mut plain_of: HashMap<&str, Box<str>> = HashMap::new();
    batch
        .iter()
        .map(|b| {
            if let Some(p) = plain_of.get(b.as_ref()) {
                return Ok(p.clone());
            }
            let p = plain_batch(b)?;
            plain_of.insert(b, p.clone());
            Ok(p)
        })
        .collect()
}

/// The column ranges of `n` columns, `block` at a time.
fn blocks(n: usize, block: usize) -> Vec<(usize, usize)> {
    legume_numeric::matrix::utils::generate_minibatch_intervals(n, 0, Some(block.max(1)))
}

/// The scored rows (all when `None`) and their names indexed: only they
/// can take a gene, by its exact name first, else as genes match anywhere
/// (a query named `TSPAN6` against a fit's `ENSG…_TSPAN6/count/spliced`).
struct RowLookup<'a> {
    n_rows: usize,
    rows: Vec<usize>,
    exact: HashMap<&'a str, usize>,
    index: GeneIndex,
}

impl<'a> RowLookup<'a> {
    fn new(names: &'a [Box<str>], scored: Option<&[usize]>) -> Self {
        let rows: Vec<usize> = scored.map_or_else(|| (0..names.len()).collect(), <[usize]>::to_vec);
        Self {
            n_rows: names.len(),
            exact: rows.iter().map(|&r| (names[r].as_ref(), r)).collect(),
            index: GeneIndex::build(&rows.iter().map(|&r| names[r].clone()).collect::<Vec<_>>()),
            rows,
        }
    }

    /// The row a functional's gene is on; `None` counts it unexpressed.
    fn row_of(&self, g: &str) -> Option<usize> {
        self.exact
            .get(g)
            .copied()
            .or_else(|| self.index.match_gene(g).map(|i| self.rows[i]))
    }
}

/// Each row's mean log-normalised expression over the cells (at most
/// about [`MEAN_CELLS`] of them, whole blocks evenly spread).
fn gene_means(
    data: &dyn ColumnBlocks,
    n_rows: usize,
    scored: Option<&[usize]>,
    block: usize,
) -> anyhow::Result<Vec<f64>> {
    // The means only rank genes into the controls' expression bins, which
    // a few tens of thousands of cells settle: every `stride`-th block.
    // One running sum per worker, merged pairwise: memory is a few of them,
    // not one per block.
    let all_blocks = blocks(data.n_cols(), block);
    let stride = all_blocks
        .len()
        .div_ceil(MEAN_CELLS.div_ceil(block.max(1)))
        .max(1);
    let all = all_blocks
        .into_par_iter()
        .step_by(stride)
        .try_fold(
            || GeneMeans::new(n_rows, scored),
            |mut m, (lb, ub)| -> anyhow::Result<GeneMeans> {
                m.add_block(&data.block(lb, ub)?);
                Ok(m)
            },
        )
        .try_reduce(
            || GeneMeans::new(n_rows, scored),
            |mut a, b| {
                a.merge(&b);
                Ok(a)
            },
        )?;
    Ok(all.means())
}

/// Every cell's score on each functional, cells × functionals. A
/// functional's gene is found by its exact name, else as genes match
/// anywhere (a query named `TSPAN6` against a fit's
/// `ENSG…_TSPAN6/count/spliced`); one found nowhere counts as unexpressed.
fn score(
    data: &dyn ColumnBlocks,
    functionals: &[Functional],
    lookup: &RowLookup<'_>,
    scored: Option<&[usize]>,
    block: usize,
) -> anyhow::Result<DMatrix<f32>> {
    let n = data.n_cols();
    anyhow::ensure!(n > 0, "--covariate: no cells to score");
    let weights = ScoreWeights::new(functionals, lookup.n_rows, |g| lookup.row_of(g), scored);
    let parts = blocks(n, block)
        .into_par_iter()
        .map(|(lb, ub)| Ok((lb, weights.score_block(&data.block(lb, ub)?))))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let mut out = DMatrix::<f32>::zeros(n, functionals.len());
    for (lb, part) in parts {
        out.rows_mut(lb, part.nrows()).copy_from(&part);
    }
    Ok(out)
}
