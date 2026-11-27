//! Covariate gene sets: per-cell programmes (the cell cycle, mitochondrial
//! and ribosomal shares, sex, dissociation stress, hypoxia, interferon)
//! that a single-cell fit can score and model.
//!
//! A covariate is one or more named gene sets, each scored one way
//! ([`SetKind`]): its share of a cell's counts, or a module score (mean log
//! expression minus that of expression-matched controls). The built-in
//! table, `data/covariate_gene_sets.tsv`, cites each source in its header;
//! the cell cycle comes from [`crate::cell_cycle`] (Tirosh et al. 2016) as
//! the covariate `cycle` with sets `S` and `G2M`.
//!
//! A table of one's own in the same form (tab-separated `covariate`, `set`,
//! `kind`, `symbol` and optionally `current`, `#` lines comments) adds
//! covariates or replaces built-in ones of the same name
//! ([`CovariateSets::with`]). Symbols match as
//! [`crate::name_matching::GeneIndex`] matches: case-insensitively, so the
//! human lists also find mouse orthologues that share a symbol, and the
//! mouse dissociation list finds human genes.

use crate::cell_cycle::{CycleGenes, CyclePhase};
use crate::name_matching::GeneIndex;
use std::path::Path;

/// The bundled table of covariate gene sets (sources in its header).
pub const COVARIATE_GENE_SETS_TSV: &str = include_str!("../data/covariate_gene_sets.tsv");

/// How a gene set is scored in a cell: `Share`, the set's share of the
/// cell's counts, or `Module`, mean log expression of the set minus that of
/// expression-matched control genes (legume-numeric's scores).
pub use legume_numeric::matrix::covariate::ScoreKind as SetKind;

/// `share` or `module`, any case.
fn parse_kind(s: &str) -> Option<SetKind> {
    match s.to_ascii_lowercase().as_str() {
        "share" => Some(SetKind::Share),
        "module" => Some(SetKind::Module),
        _ => None,
    }
}

/// The rows of `index` the genes `(symbol, current)` match, each once, in
/// gene order: a gene by its listed symbol, else its current one.
pub(crate) fn matched_symbols<'a>(
    index: &GeneIndex,
    genes: impl Iterator<Item = (&'a str, Option<&'a str>)>,
) -> Vec<usize> {
    let mut rows: Vec<usize> = Vec::new();
    for (symbol, current) in genes {
        let hit = index
            .match_gene(symbol)
            .or_else(|| current.and_then(|c| index.match_gene(c)));
        if let Some(r) = hit.filter(|r| !rows.contains(r)) {
            rows.push(r);
        }
    }
    rows
}

/// A gene of a set: its symbol as listed, and its current HGNC symbol when
/// renamed; either matches.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SetGene {
    pub symbol: Box<str>,
    pub current: Option<Box<str>>,
}

/// One scored gene set of a covariate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneSet {
    pub covariate: Box<str>,
    pub name: Box<str>,
    pub kind: SetKind,
    pub genes: Vec<SetGene>,
}

impl GeneSet {
    /// `covariate/set`, or the covariate's name alone when they agree.
    #[must_use]
    pub fn label(&self) -> String {
        if self.name == self.covariate {
            self.covariate.to_string()
        } else {
            format!("{}/{}", self.covariate, self.name)
        }
    }

    /// The rows of `index` its genes match, each once, in listed order: indices
    /// into the names `index` was built from, so a caller that built it over a
    /// subset of its rows maps them back itself.
    #[must_use]
    pub fn matched_rows(&self, index: &GeneIndex) -> Vec<usize> {
        matched_symbols(
            index,
            self.genes
                .iter()
                .map(|g| (g.symbol.as_ref(), g.current.as_deref())),
        )
    }
}

/// Covariates, each a few gene sets, in table order.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CovariateSets {
    pub sets: Vec<GeneSet>,
    pub source: Box<str>,
}

impl CovariateSets {
    /// The built-in covariates: `cycle` (Tirosh et al. 2016) and the
    /// bundled table.
    ///
    /// # Panics
    /// Never: the bundled tables are checked by this crate's tests.
    #[must_use]
    pub fn builtin() -> Self {
        let mut sets = cycle_sets(&CycleGenes::tirosh_2016());
        let table = Self::parse(COVARIATE_GENE_SETS_TSV, "covariate_gene_sets.tsv")
            .expect("the bundled table parses");
        sets.extend(table.sets);
        Self {
            sets,
            source: "built-in (Tirosh 2016; covariate_gene_sets.tsv)".into(),
        }
    }

    /// Parse a table: tab-separated `covariate`, `set`, `kind`, `symbol`
    /// and optionally `current`; blank and `#` lines and a header line
    /// (first column `covariate`) skipped.
    ///
    /// # Errors
    /// A line missing a column or with a kind that is neither `share` nor
    /// `module`; one set listed with two kinds; or no genes at all.
    pub fn parse(text: &str, source: &str) -> Result<Self, String> {
        let mut sets: Vec<GeneSet> = Vec::new();
        for (k, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').map(str::trim).collect();
            if cols[0].eq_ignore_ascii_case("covariate") {
                continue;
            }
            let at = |what: &str| format!("{source} line {}: {what}", k + 1);
            let [covariate, set, kind, symbol] =
                [0, 1, 2, 3].map(|i| cols.get(i).copied().unwrap_or(""));
            if covariate.is_empty() || set.is_empty() || symbol.is_empty() {
                return Err(at("needs covariate, set, kind and symbol"));
            }
            let kind = parse_kind(kind)
                .ok_or_else(|| at(&format!("kind {kind:?} is not share or module")))?;
            let current = cols
                .get(4)
                .filter(|c| !c.is_empty() && !c.eq_ignore_ascii_case(symbol))
                .map(|c| Box::from(*c));
            let gene = SetGene {
                symbol: symbol.into(),
                current,
            };
            match sets
                .iter_mut()
                .find(|s| *s.covariate == *covariate && *s.name == *set)
            {
                Some(s) if s.kind != kind => {
                    return Err(at(&format!("set {covariate}/{set} listed as two kinds")))
                }
                Some(s) => s.genes.push(gene),
                None => sets.push(GeneSet {
                    covariate: covariate.into(),
                    name: set.into(),
                    kind,
                    genes: vec![gene],
                }),
            }
        }
        if sets.is_empty() {
            return Err(format!("{source}: no gene sets"));
        }
        Ok(Self {
            sets,
            source: source.into(),
        })
    }

    /// [`CovariateSets::parse`] of the file at `path`.
    ///
    /// # Errors
    /// A file that cannot be read, or as [`CovariateSets::parse`].
    pub fn read(path: &Path) -> Result<Self, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        Self::parse(&text, &path.display().to_string())
    }

    /// These covariates with `other`'s: a covariate `other` names replaces
    /// the one here of that name (all its sets), and a new one is added.
    #[must_use]
    pub fn with(self, other: Self) -> Self {
        let replaced: Vec<&str> = other.covariates();
        let mut sets: Vec<GeneSet> = self
            .sets
            .iter()
            .filter(|s| !replaced.contains(&s.covariate.as_ref()))
            .cloned()
            .collect();
        sets.extend(other.sets.iter().cloned());
        Self {
            sets,
            source: format!("{} with {}", self.source, other.source).into(),
        }
    }

    /// The covariates' names, each once, in table order.
    #[must_use]
    pub fn covariates(&self) -> Vec<&str> {
        let mut names: Vec<&str> = Vec::new();
        for s in &self.sets {
            if !names.contains(&s.covariate.as_ref()) {
                names.push(&s.covariate);
            }
        }
        names
    }

    /// The sets of `covariate`.
    pub fn of<'a>(&'a self, covariate: &'a str) -> impl Iterator<Item = &'a GeneSet> + 'a {
        self.sets.iter().filter(move |s| *s.covariate == *covariate)
    }
}

/// The cell cycle as covariate `cycle`: module sets `S` and `G2M`.
fn cycle_sets(cycle: &CycleGenes) -> Vec<GeneSet> {
    [(CyclePhase::S, "S"), (CyclePhase::G2M, "G2M")]
        .into_iter()
        .map(|(phase, name)| GeneSet {
            covariate: "cycle".into(),
            name: name.into(),
            kind: SetKind::Module,
            genes: cycle
                .genes
                .iter()
                .filter(|g| g.phase == phase)
                .map(|g| SetGene {
                    symbol: g.symbol.clone(),
                    current: g.current.clone(),
                })
                .collect(),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_builtin_table_has_every_covariate_with_its_sets_and_kinds() {
        let b = CovariateSets::builtin();
        assert_eq!(
            b.covariates(),
            [
                "cycle",
                "mito",
                "ribo",
                "sex",
                "stress",
                "hypoxia",
                "interferon"
            ]
        );
        let count = |c: &str, s: &str| b.of(c).find(|x| &*x.name == s).map_or(0, |x| x.genes.len());
        assert_eq!(count("cycle", "S"), 43);
        assert_eq!(count("cycle", "G2M"), 54);
        assert_eq!(count("mito", "mito"), 37);
        assert_eq!(count("ribo", "ribo"), 87);
        assert_eq!((count("sex", "female"), count("sex", "male")), (1, 4));
        assert_eq!(count("stress", "dissociation"), 140);
        assert_eq!(count("hypoxia", "hypoxia"), 200);
        assert_eq!(
            (count("interferon", "alpha"), count("interferon", "gamma")),
            (97, 200)
        );
        assert!(b.of("mito").all(|s| s.kind == SetKind::Share));
        assert!(b.of("ribo").all(|s| s.kind == SetKind::Share));
        assert!(b.of("hypoxia").all(|s| s.kind == SetKind::Module));
        // RPS4Y1 marks sex, not ribosomes.
        assert!(!b
            .of("ribo")
            .any(|s| s.genes.iter().any(|g| &*g.symbol == "RPS4Y1")));
    }

    #[test]
    fn a_table_of_ones_own_replaces_or_adds_covariates() {
        let mine = CovariateSets::parse(
            "covariate\tset\tkind\tsymbol\nstress\theat\tmodule\tHSPA1A\nstress\theat\tmodule\tHSPA1B\nown\tx\tshare\tGAPDH\n",
            "mine",
        )
        .unwrap();
        let b = CovariateSets::builtin().with(mine);
        let stress: Vec<_> = b.of("stress").collect();
        assert_eq!(stress.len(), 1);
        assert_eq!(&*stress[0].name, "heat");
        assert_eq!(stress[0].genes.len(), 2);
        assert_eq!(b.covariates().last(), Some(&"own"));
        assert_eq!(b.of("own").next().unwrap().label(), "own/x");
        assert_eq!(b.of("mito").next().unwrap().label(), "mito");
    }

    #[test]
    fn a_bad_table_says_which_line() {
        assert!(CovariateSets::parse("a\tb\tweird\tG\n", "t")
            .unwrap_err()
            .contains("line 1"));
        assert!(CovariateSets::parse("a\tb\tshare\n", "t").is_err());
        assert!(
            CovariateSets::parse("a\tb\tshare\tG\na\tb\tmodule\tH\n", "t")
                .unwrap_err()
                .contains("two kinds")
        );
        assert!(CovariateSets::parse("# nothing\n", "t").is_err());
    }

    #[test]
    fn sets_match_through_the_lenient_index() {
        let index = GeneIndex::build(&[
            Box::from("ENSG00000229807_XIST"),
            Box::from("Fos"),
            Box::from("ENSG00000129824_RPS4Y1/count/spliced"),
        ]);
        let b = CovariateSets::builtin();
        let female = b.of("sex").find(|s| &*s.name == "female").unwrap();
        assert_eq!(female.matched_rows(&index), [0]);
        let male = b.of("sex").find(|s| &*s.name == "male").unwrap();
        assert_eq!(male.matched_rows(&index), [2]);
        let stress = b.of("stress").next().unwrap();
        assert_eq!(stress.matched_rows(&index), [1]);
    }
}
