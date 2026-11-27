//! Cell-cycle marker genes: auxiliary gene sets for scoring each cell's
//! S and G2/M programmes.
//!
//! The default sets, bundled as `data/cell_cycle_tirosh2016.tsv`, are the
//! 43 S-phase and 54 G2/M-phase genes of
//!
//! > Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular
//! > ecosystem of metastatic melanoma by single-cell RNA-seq.
//! > *Science* 352(6282):189–196 (2016). doi:10.1126/science.aad0501
//!
//! in the form Seurat distributes them (`cc.genes`), each with its current
//! HGNC symbol where it has been renamed (Seurat's `cc.genes.updated.2019`:
//! MLF1IP → CENPU, FAM64A → PIMREG, HN1 → JPT1); either spelling matches.
//!
//! The lists are human HGNC symbols. Matched case-insensitively (as
//! [`crate::name_matching::GeneIndex`] matches), they also find the mouse
//! orthologues these conserved genes share a symbol with (`Mki67`,
//! `Top2a`, `Pcna`); a dictionary that names them otherwise, or a user who
//! wants other genes, reads a table of its own in the same form
//! ([`CycleGenes::parse`], [`CycleGenes::read`]): tab-separated `phase`
//! (`S` or `G2M`), `symbol` and optionally `current`, `#` lines comments.

use crate::name_matching::GeneIndex;
use std::path::Path;

/// The bundled default table (Tirosh et al. 2016).
pub const TIROSH_2016_TSV: &str = include_str!("../data/cell_cycle_tirosh2016.tsv");

/// A phase of the cycle a gene set marks.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum CyclePhase {
    /// DNA synthesis.
    S,
    /// G2 and mitosis.
    G2M,
}

impl CyclePhase {
    /// The phase a table names: `S`, or `G2M` (`G2/M` too), any case.
    #[must_use]
    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_uppercase().replace('/', "").as_str() {
            "S" => Some(Self::S),
            "G2M" => Some(Self::G2M),
            _ => None,
        }
    }
}

/// One gene of a set: its symbol, and another spelling (a later HGNC
/// symbol) that names the same gene.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CycleGene {
    pub phase: CyclePhase,
    pub symbol: Box<str>,
    pub current: Option<Box<str>>,
}

/// S and G2/M gene sets: the bundled default or a user's own.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CycleGenes {
    pub genes: Vec<CycleGene>,
    /// Where they came from, for a log line or a table's metadata.
    pub source: Box<str>,
}

impl CycleGenes {
    /// The default sets: Tirosh et al. 2016.
    ///
    /// # Panics
    /// Never: the bundled table is checked by the tests.
    #[must_use]
    pub fn tirosh_2016() -> Self {
        Self::parse(TIROSH_2016_TSV, "Tirosh et al. 2016 (Science 352:189)")
            .expect("the bundled cell-cycle table parses")
    }

    /// A table in the bundled form: tab-separated `phase`, `symbol` and an
    /// optional `current`; `#` lines and a `phase` header skipped.
    ///
    /// # Errors
    /// A line with no symbol or a phase that is neither `S` nor `G2M`, or
    /// a table missing either phase.
    pub fn parse(text: &str, source: &str) -> Result<Self, String> {
        let mut genes = Vec::new();
        for (k, line) in text.lines().enumerate() {
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            let mut cols = line.split('\t').map(str::trim);
            let (phase, symbol) = (cols.next().unwrap_or(""), cols.next().unwrap_or(""));
            if phase.eq_ignore_ascii_case("phase") {
                continue;
            }
            let phase = CyclePhase::parse(phase).ok_or_else(|| {
                format!("{source} line {}: phase {phase:?} is not S or G2M", k + 1)
            })?;
            if symbol.is_empty() {
                return Err(format!("{source} line {}: no symbol", k + 1));
            }
            let current = cols
                .next()
                .filter(|c| !c.is_empty() && !c.eq_ignore_ascii_case(symbol))
                .map(Box::from);
            genes.push(CycleGene {
                phase,
                symbol: symbol.into(),
                current,
            });
        }
        for p in [CyclePhase::S, CyclePhase::G2M] {
            if !genes.iter().any(|g| g.phase == p) {
                return Err(format!("{source}: no {p:?} genes"));
            }
        }
        Ok(Self {
            genes,
            source: source.into(),
        })
    }

    /// [`CycleGenes::parse`] of the file at `path`.
    ///
    /// # Errors
    /// A file that cannot be read, or as [`CycleGenes::parse`].
    pub fn read(path: &Path) -> Result<Self, String> {
        let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
        Self::parse(&text, &path.display().to_string())
    }

    /// The symbols of `phase`, as listed.
    pub fn symbols(&self, phase: CyclePhase) -> impl Iterator<Item = &str> {
        self.genes
            .iter()
            .filter(move |g| g.phase == phase)
            .map(|g| &*g.symbol)
    }

    /// The rows of `index` that `phase`'s genes match, by either spelling,
    /// each row once, in the list's order.
    #[must_use]
    pub fn matched_rows(&self, index: &GeneIndex, phase: CyclePhase) -> Vec<usize> {
        crate::gene_sets::matched_symbols(
            index,
            self.genes
                .iter()
                .filter(|g| g.phase == phase)
                .map(|g| (g.symbol.as_ref(), g.current.as_deref())),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_default_is_tirosh_s_43_and_54_genes() {
        let c = CycleGenes::tirosh_2016();
        assert_eq!(c.symbols(CyclePhase::S).count(), 43);
        assert_eq!(c.symbols(CyclePhase::G2M).count(), 54);
        let renamed: Vec<_> = c
            .genes
            .iter()
            .filter_map(|g| g.current.as_deref())
            .collect();
        assert_eq!(renamed, ["CENPU", "PIMREG", "JPT1"]);
        assert!(c
            .symbols(CyclePhase::S)
            .all(|g| !c.symbols(CyclePhase::G2M).any(|h| h == g)));
    }

    #[test]
    fn a_mouse_dictionary_and_a_renamed_symbol_match_once() {
        let names: Vec<Box<str>> = ["Mki67", "Top2a", "Pcna", "Actb", "CENPU", "JPT1"]
            .iter()
            .map(|s| Box::from(*s))
            .collect();
        let index = GeneIndex::build(&names);
        let c = CycleGenes::tirosh_2016();
        assert_eq!(c.matched_rows(&index, CyclePhase::S), vec![2, 4]);
        let g2m = c.matched_rows(&index, CyclePhase::G2M);
        assert!(g2m.contains(&0) && g2m.contains(&1) && g2m.contains(&5));
        assert!(!g2m.contains(&3));
    }

    #[test]
    fn lenient_names_match_but_near_names_do_not() {
        // Ensembl-joined and aux-suffixed rows match; HN1L, MKI67IP and
        // PCNA-AS1 are other genes and must not.
        let names: Vec<Box<str>> = [
            "ENSG00000148773_MKI67",
            "TOP2A/count/spliced",
            "HN1L",
            "MKI67IP",
            "PCNA-AS1",
        ]
        .iter()
        .map(|s| Box::from(*s))
        .collect();
        let index = GeneIndex::build(&names);
        let c = CycleGenes::tirosh_2016();
        let g2m = c.matched_rows(&index, CyclePhase::G2M);
        assert_eq!(g2m, vec![1, 0]);
        assert!(c.matched_rows(&index, CyclePhase::S).is_empty());
    }

    #[test]
    fn a_user_table_overrides_and_a_bad_one_is_refused() {
        let own = CycleGenes::parse("# mine\nS\tFoo\nG2/M\tBar\tBaz\n", "mine").unwrap();
        assert_eq!(own.symbols(CyclePhase::S).collect::<Vec<_>>(), ["Foo"]);
        assert_eq!(own.genes[1].current.as_deref(), Some("Baz"));
        assert!(CycleGenes::parse("S\tFoo\n", "x")
            .unwrap_err()
            .contains("no G2M"));
        assert!(CycleGenes::parse("G1\tFoo\n", "x")
            .unwrap_err()
            .contains("not S or G2M"));
    }
}
