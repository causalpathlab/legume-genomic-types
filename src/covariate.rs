//! Per-cell covariates as a second axis of batch, one implementation for
//! every tool that models batch by a per-gene δ (senna, pinto, chickpea):
//! each chosen covariate's gene sets ([`crate::gene_sets`]) scored in every
//! cell, the scores binned together (legume-numeric's
//! [`legume_numeric::matrix::covariate`]: centred within batches,
//! standardised robustly, one Gaussian mixture, k by ICL), and each cell's
//! batch label made `batch|cvK` ([`crate::batch_label`]), so the δ is fitted
//! per batch and bin with no change to the model.
//!
//! A tool supplies its counts as [`ColumnBlocks`], its cells' plain batch
//! labels and a random projection of them (for the groups of like cells
//! the one-bin guard looks at), and gets back each cell's label and the
//! fitted [`CovariateModel`], which bins new cells the same way. Where the
//! labels apply, how the projection is made and what is written beside a
//! run stay the tool's.
//!
//! - `--covariate all` (the default) bins the covariates whose genes match
//!   and that show groups of their own; a list bins exactly those; `none`
//!   models nothing. A run recorded before these flags replays as `none`
//!   ([`CovariateArgs`]'s `Default`).

mod fit;
mod outputs;

pub use fit::{model_covariates, CovariateModel, Modelled};
pub use outputs::{write_cells, write_covariate_effect};

use crate::batch_label;
use clap::Args;
use legume_numeric::matrix::gmm::Criterion;
use nalgebra_sparse::CscMatrix;
use serde::{Deserialize, Serialize};

/// The label part of a covariate bin: `cv0`, `cv1`, ….
pub const AXIS: &str = "cv";

/// Where a run's covariate outputs sit beside it (`prefix` its output
/// prefix), named alike in every tool.
pub mod files {
    /// The fitted [`super::CovariateModel`], JSON.
    #[must_use]
    pub fn model(prefix: &str) -> String {
        format!("{prefix}.covariate_model.json")
    }

    /// Every cell's scores and label ([`super::write_cells`]).
    #[must_use]
    pub fn cells(prefix: &str) -> String {
        format!("{prefix}.covariates.parquet")
    }

    /// Each gene's covariate effect ([`super::write_covariate_effect`]).
    #[must_use]
    pub fn effect(prefix: &str) -> String {
        format!("{prefix}.covariate_effect.parquet")
    }

    /// The `i`-th data file's batch labels, one line per column.
    #[must_use]
    pub fn batch(prefix: &str, i: usize) -> String {
        format!("{prefix}.covariate.batch.{i}.txt")
    }
}

/// A tool's counts, rows (features) × columns (cells), read a block of
/// columns at a time. Every row is returned; the covariate module restricts
/// to the rows it scores.
pub trait ColumnBlocks: Sync {
    fn n_cols(&self) -> usize;
    fn row_names(&self) -> anyhow::Result<Vec<Box<str>>>;
    /// Columns `lb..ub`, all rows.
    fn block(&self, lb: usize, ub: usize) -> anyhow::Result<CscMatrix<f32>>;
}

/// The flags a fit takes for its covariates.
#[derive(Args, Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct CovariateArgs {
    #[arg(
        long = "covariate",
        value_name = "NAMES",
        value_delimiter = ',',
        default_value = "all",
        help = "Per-cell covariates binned as a second axis of batch: all, none, or a list \
                (cycle, mito, ribo, sex, stress, hypoxia, interferon)",
        long_help = "Per-cell covariates binned as a second axis of batch (labels batch|cvK).\n\
                     all (default): every built-in covariate whose genes match the features,\n\
                     binning those with groups of their own (each judged on its scores\n\
                     together); none: nothing; or a list, binned as given:\n\
                     cycle (S, G2/M; Tirosh 2016), mito, ribo (count shares; HGNC), sex\n\
                     (XIST against Y genes; Staedtler 2013), stress (dissociation genes;\n\
                     van den Brink 2017), hypoxia, interferon (alpha, gamma; MSigDB\n\
                     Hallmark), or a covariate named in --covariate-genes."
    )]
    pub covariates: Vec<Box<str>>,

    #[arg(
        long,
        value_name = "FILE",
        help = "Gene sets of your own (covariate, set, kind share|module, symbol): add \
                covariates or replace built-in ones of the same name"
    )]
    pub covariate_genes: Option<Box<str>>,

    #[arg(
        long,
        help = "With --covariate: this many bins [default: chosen by --covariate-criterion]"
    )]
    pub covariate_bins: Option<usize>,

    #[arg(
        long,
        default_value_t = 10,
        help = "With --covariate: most bins the criterion may choose"
    )]
    pub covariate_bins_max: usize,

    #[arg(
        long,
        value_enum,
        default_value_t = Criterion::Icl,
        help = "With --covariate: what chooses the number of bins: icl (BIC plus the \
                memberships' entropy, which charges overlapping bins) or bic"
    )]
    pub covariate_criterion: Criterion,
}

/// A run recorded before covariates existed replays without them.
impl Default for CovariateArgs {
    fn default() -> Self {
        Self {
            covariates: vec!["none".into()],
            covariate_genes: None,
            covariate_bins: None,
            covariate_bins_max: 10,
            covariate_criterion: Criterion::Icl,
        }
    }
}

/// What `--covariate` asked for.
pub(crate) enum Selection {
    None,
    All,
    Named(Vec<Box<str>>),
}

impl CovariateArgs {
    pub(crate) fn selection(&self) -> anyhow::Result<Selection> {
        let words: Vec<&str> = self.covariates.iter().map(|w| w.trim()).collect();
        match words.as_slice() {
            [] | ["none"] => Ok(Selection::None),
            ["all"] => Ok(Selection::All),
            many => {
                anyhow::ensure!(
                    !many.iter().any(|w| *w == "all" || *w == "none"),
                    "--covariate: all or none stands alone, not in a list"
                );
                Ok(Selection::Named(many.iter().map(|w| (*w).into()).collect()))
            }
        }
    }

    /// Whether covariates are to be modelled, `false` for `none`.
    ///
    /// # Errors
    /// A malformed `--covariate`.
    pub fn is_on(&self) -> anyhow::Result<bool> {
        Ok(!matches!(self.selection()?, Selection::None))
    }

    /// Where covariates do not yet work (`what`): the default `all` says so
    /// in the log and models none; a list asked for by name is refused.
    ///
    /// # Errors
    /// Covariates named, or a malformed `--covariate`.
    pub fn unsupported(&self, what: &str) -> anyhow::Result<()> {
        self.enabled(&[(true, what)]).map(|_| ())
    }

    /// Whether covariates are modelled in a load that may not support them:
    /// each `(applies, what)` that applies keeps `all` off (said in the log)
    /// and refuses a list; `false` too for `none`.
    ///
    /// # Errors
    /// Covariates named where one applies, or a malformed `--covariate`.
    pub fn enabled(&self, blockers: &[(bool, &str)]) -> anyhow::Result<bool> {
        let selection = self.selection()?;
        let mut on = !matches!(selection, Selection::None);
        for &(_, what) in blockers.iter().filter(|(applies, _)| *applies) {
            match selection {
                Selection::None => {}
                Selection::All => log::info!("--covariate all: not modelled with {what}"),
                Selection::Named(_) => {
                    anyhow::bail!("--covariate does not yet work with {what}")
                }
            }
            on = false;
        }
        Ok(on)
    }
}

/// `label` without the covariate axis an earlier fit gave it: `rep1|cv0` is
/// `rep1`, `rep1` stays.
///
/// # Errors
/// A label with any other second axis, which a covariate fit cannot extend.
pub fn plain_batch(label: &str) -> anyhow::Result<Box<str>> {
    match batch_label::split(label).as_slice() {
        [b] => Ok((*b).into()),
        [b, c] if bin_of(c).is_some() => Ok((*b).into()),
        _ => anyhow::bail!(
            "--covariate: batch label {label:?} already has more than one part (joined by '{}'); \
             give plain batch labels or --covariate none",
            batch_label::AXIS_SEP
        ),
    }
}

/// The bin a label part names, `cv0` → 0; `None` for any other part.
#[must_use]
pub fn bin_of(part: &str) -> Option<usize> {
    part.strip_prefix(AXIS)
        .filter(|k| !k.is_empty() && k.bytes().all(|b| b.is_ascii_digit()))?
        .parse()
        .ok()
}

/// `batch|cvK`, or the plain `batch` without a bin.
///
/// # Errors
/// A batch label that cannot be a label part.
pub fn compose(batch: &str, bin: Option<usize>) -> anyhow::Result<Box<str>> {
    Ok(match bin {
        Some(c) => batch_label::compose(&[batch, &format!("{AXIS}{c}")])
            .map_err(anyhow::Error::msg)?
            .into(),
        None => batch.into(),
    })
}

#[cfg(test)]
mod tests;
