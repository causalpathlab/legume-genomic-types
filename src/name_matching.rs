//! Matching gene names across naming conventions: a [`GeneIndex`] over a
//! dictionary's row names resolves HGNC symbols, Ensembl ids, `ENSG…_SYMBOL`
//! joins, renamed HGNC symbols and genomic loci to rows.
//!
//! Moved from data-beans' `utilities::name_matching` (which re-exports it)
//! so crates below data-beans can match gene names without depending on it.

use crate::coordinates::{is_locus, locus_key};
use rayon::prelude::*;
use rustc_hash::FxHashMap as HashMap;

/// The separator joining a feature id and its name in a row name
/// (`ENSG00000153563_CD8A`).
pub const ROW_SEP: &str = "_";

/// Below this matched-gene fraction, a gene-set reconciliation warns rather
/// than logs at info: a thin overlap silently produces no enrichment.
pub const COVERAGE_WARN_FRAC: f32 = 0.5;

/// Heuristic: a lower-cased Ensembl-style stable id (`ensg…`, `ensmusg…`,
/// `enst…`). Used to index/look up the *leading* id segment of an
/// `ENSG…_SYMBOL` name so bare-`ENSG` and `ENSG_SYMBOL` forms reconcile both
/// ways without a linear scan.
fn is_ensembl_id(s: &str) -> bool {
    s.len() >= 8 && s.starts_with("ens") && s.bytes().any(|b| b.is_ascii_digit())
}

/// Curated HGNC **old symbol → current symbol** renames, lower-cased.
///
/// A symbol match is exact, so a marker panel written against an older HGNC release
/// silently loses every gene HGNC has since renamed — the gene is in the matrix under its
/// new name, but the panel asks for the old one and gets nothing back. That is invisible in
/// the output: the type just scores on fewer genes (or is dropped entirely). This table is
/// what closes the gap.
///
/// It is **curated, not exhaustive** — the families that actually recur in single-cell
/// marker panels (histones, the `MARCH`/`SEPT` families that Excel also mangles into dates,
/// the selenoproteins, and the well-known one-off renames). Systematic families are handled
/// by rule in [`alias_candidates`] rather than enumerated here. Entries are one-directional
/// in this table but matched **both ways** at lookup, so it does not matter whether the
/// matrix or the panel is the one carrying the old name.
static HGNC_RENAMES: &[(&str, &str)] = &[
    // Histones — HGNC's 2019 systematic renaming; heavily used as cell-cycle / S-phase
    // markers, so a stale panel loses much of its S-phase signature.
    ("h1f0", "h1-0"),
    ("h1fx", "h1-10"),
    ("hist1h1b", "h1-5"),
    ("hist1h1c", "h1-2"),
    ("hist1h1d", "h1-3"),
    ("hist1h1e", "h1-4"),
    ("hist1h2ac", "h2ac6"),
    ("hist1h2bk", "h2bc12"),
    ("hist1h4c", "h4c3"),
    ("hist2h2be", "h2bc21"),
    ("hist3h2a", "h2ac25"),
    ("h2afx", "h2ax"),
    ("h2afv", "h2az2"),
    ("h2afz", "h2az1"),
    ("h2afy", "macroh2a1"),
    ("h3f3a", "h3-3a"),
    ("h3f3b", "h3-3b"),
    // Mitochondrial amidoxime-reducing components (note: NOT the MARCH family below).
    ("marc1", "mtarc1"),
    ("marc2", "mtarc2"),
    // Selenoproteins.
    ("sepp1", "selenop"),
    ("selt", "selenot"),
    ("sepw1", "selenow"),
    // One-off renames common in immune / proliferation panels.
    ("fam129a", "niban1"),
    ("fam129b", "niban2"),
    ("fam129c", "niban3"),
    ("rarres3", "plaat4"),
    ("fyb", "fyb1"),
    ("cd97", "adgre5"),
    ("gpr56", "adgrg1"),
    ("kiaa0101", "pclaf"),
    ("c10orf54", "vsir"),
    ("tmem66", "saraf"),
    ("atpif1", "atp5if1"),
    ("fam46c", "tent5c"),
    ("whsc1", "nsd2"),
];

/// `HGNC_RENAMES` as a lookup: old → new (`fwd`) and new → old (`rev`).
///
/// Kept as two maps rather than one seeded in both directions. The table's key sets happen to
/// be disjoint today, so one map would give identical answers — but the moment someone adds a
/// *chained* rename (`A→B` alongside an existing `B→C`), a single map has two entries for `B`
/// and silently keeps whichever was inserted last. Two maps cannot lose that way.
struct RenameMaps {
    fwd: HashMap<&'static str, &'static str>,
    rev: HashMap<&'static str, &'static str>,
}

/// The lazily-built rename lookups. Built once; every `match_gene` miss consults them.
fn rename_maps() -> &'static RenameMaps {
    static MAPS: std::sync::OnceLock<RenameMaps> = std::sync::OnceLock::new();
    MAPS.get_or_init(|| RenameMaps {
        fwd: HGNC_RENAMES.iter().copied().collect(),
        rev: HGNC_RENAMES.iter().map(|&(o, n)| (n, o)).collect(),
    })
}

/// The numeric suffix of a `{prefix}{n}` symbol (`numeric_suffix("march12", "march") == 12`),
/// or `None` if `sym` does not have exactly that shape.
fn numeric_suffix(sym: &str, prefix: &str) -> Option<u32> {
    sym.strip_prefix(prefix)
        .filter(|d| !d.is_empty() && d.bytes().all(|b| b.is_ascii_digit()))
        .and_then(|d| d.parse().ok())
}

/// Alternative HGNC symbols for `sym` (already lower-cased): the table above in both
/// directions, plus the two rule-based families whose members are too numerous to enumerate
/// and whose rename is purely mechanical — `MARCH{n}` ↔ `MARCHF{n}` (membrane-associated
/// ring-CH E3 ligases) and `SEPT{n}` ↔ `SEPTIN{n}` (septins). Both families were renamed
/// precisely because spreadsheets kept coercing them to dates, so panels in the wild carry
/// either form.
///
/// `MARCH1`/`MARC1` do not collide: `MARC1` is in the table (→ `MTARC1`) and the `march`
/// rule only fires on the literal `march` prefix.
fn alias_candidates(sym: &str) -> Vec<String> {
    let maps = rename_maps();
    let mut out = Vec::new();
    if let Some(&new) = maps.fwd.get(sym) {
        out.push(new.to_string());
    }
    if let Some(&old) = maps.rev.get(sym) {
        out.push(old.to_string());
    }
    for (old, new) in [("march", "marchf"), ("sept", "septin")] {
        if let Some(n) = numeric_suffix(sym, old) {
            out.push(format!("{new}{n}"));
        }
        if let Some(n) = numeric_suffix(sym, new) {
            out.push(format!("{old}{n}"));
        }
    }
    out
}

/// Pre-built index over a gene-name vocabulary for fast marker→row matching.
/// Resolves a query gene in tiers, returning the first matching row:
///   0. a locus query (`chr:start-end`) matches only a locus row, by its
///      locus key with the chromosome case kept, and stops here: loci never
///      go through the case-insensitive or fuzzy tiers below,
///   1. exact (case-insensitive) full-name match,
///   2. last `_`-segment symbol match (`CD8A` ↔ `ENSG…_CD8A`),
///   3. leading Ensembl-id segment match (`ENSG…` ↔ `ENSG…_CD8A`),
///   4. decompose a combined `ENSG…_SYMBOL` query and retry tiers 2–3 per part,
///   5. HGNC alias retry ([`alias_candidates`]: `HIST1H4C` ↔ `H4C3`, `MARCH2` ↔ `MARCHF2`),
///   6. fallback to a prefix / `_x_` segment match.
///
/// Tiers 1–5 are O(1) hash lookups, and the fallback tests only the rows that
/// hold the query's first `_`-segment. Build once, match many — replaces the O(genes × markers)
/// `.position(..)` scan over every name. Tier 1 preferring an exact match
/// over an earlier-indexed suffix match is the one intended refinement vs a
/// pure positional scan. Tiers 3–4 make HGNC / ENSG / `ENSG_HGNC` reconcile
/// in either direction (gene-set sources mix these conventions), and tier 5
/// does the same across HGNC *releases* — the matrix and the marker panel are
/// routinely built against different ones.
pub struct GeneIndex {
    lowered: Vec<String>,
    exact: HashMap<String, usize>,
    symbol: HashMap<String, usize>,
    ensg: HashMap<String, usize>,
    /// Every `_`-separated segment of a lowered name, with the rows that hold
    /// it in ascending order: the candidates of the fallback tier.
    segment: HashMap<String, Vec<usize>>,
    /// Rows with a locus part ([`locus_part`]), by their name as written,
    /// then by locus key plus the rest of the name, then (for rows with a
    /// rest) by the bare locus key. Loci match only here, case kept.
    locus_raw: HashMap<Box<str>, usize>,
    locus: HashMap<Box<str>, usize>,
}

/// A name's locus part and the rest after it: the whole name when it is a
/// locus, its `/`-core (`chr1:1-2/count/spliced`), or the locus before an
/// `id{ROW_SEP}name` join (`chr1:1-2_GENE1`). `None` when no part is a locus.
pub fn locus_part(name: &str) -> Option<(&str, &str)> {
    if is_locus(name) {
        return Some((name, ""));
    }
    let core = name.split('/').next().unwrap_or(name);
    let cut = if is_locus(core) {
        core.len()
    } else {
        let colon = core.find(':')?;
        colon + 1 + core[colon + 1..].find(ROW_SEP)?
    };
    let (part, rest) = name.split_at(cut);
    is_locus(part).then_some((part, rest))
}

/// The key a name with a locus part is matched by: the locus key, then the
/// rest of the name as written.
fn locus_match_key(part: &str, rest: &str) -> Option<Box<str>> {
    let key = locus_key(part)?;
    Some(if rest.is_empty() {
        key
    } else {
        format!("{key}{rest}").into_boxed_str()
    })
}

impl GeneIndex {
    /// Build the index from the dictionary's gene-name order. The first row
    /// wins on duplicate keys (matching positional-scan semantics).
    #[must_use]
    pub fn build(gene_names: &[Box<str>]) -> Self {
        // A whole-locus row gets an empty lowered name, which keeps it out
        // of every gene tier, the fallback scan included. A row with a
        // locus part plus a rest stays in the gene tiers for its full name.
        let lowered: Vec<String> = gene_names
            .par_iter()
            .map(|g| {
                if is_locus(g) {
                    String::new()
                } else {
                    g.to_lowercase()
                }
            })
            .collect();
        let mut locus_raw: HashMap<Box<str>, usize> = HashMap::default();
        let mut locus: HashMap<Box<str>, usize> = HashMap::default();
        let mut bare: Vec<(Box<str>, usize)> = Vec::new();
        for (i, g) in gene_names.iter().enumerate() {
            let Some((part, rest)) = locus_part(g) else {
                continue;
            };
            locus_raw.entry(g.clone()).or_insert(i);
            if let Some(key) = locus_match_key(part, rest) {
                locus.entry(key).or_insert(i);
            }
            if !rest.is_empty() {
                bare.extend(locus_key(part).map(|k| (k, i)));
            }
        }
        // A whole-locus row wins its key over a row that only starts with it.
        for (key, i) in bare {
            locus.entry(key).or_insert(i);
        }
        let mut exact: HashMap<String, usize> = HashMap::default();
        let mut symbol: HashMap<String, usize> = HashMap::default();
        let mut ensg: HashMap<String, usize> = HashMap::default();
        let mut segment: HashMap<String, Vec<usize>> = HashMap::default();
        for (i, low) in lowered.iter().enumerate() {
            if low.is_empty() {
                continue;
            }
            exact.entry(low.clone()).or_insert(i);
            for seg in low.split('_') {
                let rows = segment.entry(seg.to_string()).or_default();
                if rows.last() != Some(&i) {
                    rows.push(i);
                }
            }
            // Strip a faba-style aux suffix first (`SYMBOL/count/spliced` →
            // symbol is the leading `/`-segment), then an Ensembl-style prefix
            // (`ENSG…_CD8A` → symbol is the trailing `_`-segment). Handles
            // either convention, or both combined (`ENSG…_CD8A/count/spliced`).
            let core = low.split('/').next().unwrap_or(low);
            if let Some(sym) = core.rsplit('_').next() {
                symbol.entry(sym.to_string()).or_insert(i);
            }
            // Also index the *leading* segment when it is an Ensembl id, so a
            // bare `ENSG…` query resolves to an `ENSG…_SYMBOL` row (and back).
            let head = core.split('_').next().unwrap_or(core);
            if is_ensembl_id(head) {
                ensg.entry(head.to_string()).or_insert(i);
            }
        }
        Self {
            lowered,
            exact,
            symbol,
            ensg,
            segment,
            locus_raw,
            locus,
        }
    }

    /// Row index for `gene`, or `None` if unmatched (tiers above).
    #[must_use]
    pub fn match_gene(&self, gene: &str) -> Option<usize> {
        // A name with a locus part matches only a row with one: the same
        // name as written first, then by locus key plus the rest. Strictly:
        // no case folding, aliasing or prefix fallback.
        if let Some((part, rest)) = locus_part(gene) {
            if let Some(&i) = self.locus_raw.get(gene) {
                return Some(i);
            }
            return locus_match_key(part, rest).and_then(|k| self.locus.get(&k).copied());
        }
        let gl = gene.to_lowercase();
        if let Some(&i) = self.exact.get(&gl) {
            return Some(i);
        }
        if let Some(&i) = self.symbol.get(&gl) {
            return Some(i);
        }
        if let Some(&i) = self.ensg.get(&gl) {
            return Some(i);
        }
        // Decompose a combined `ENSG…_SYMBOL[/aux]` query: match its trailing
        // symbol or leading Ensembl id against the per-part indices.
        let core = gl.split('/').next().unwrap_or(&gl);
        if let Some(sym) = core.rsplit('_').next() {
            if sym != gl {
                if let Some(&i) = self.symbol.get(sym) {
                    return Some(i);
                }
            }
        }
        let head = core.split('_').next().unwrap_or(core);
        if is_ensembl_id(head) {
            if let Some(&i) = self.ensg.get(head) {
                return Some(i);
            }
        }
        // HGNC alias retry: the query and the vocabulary can be built against different HGNC
        // releases (`HIST1H4C` in the panel, `H4C3` in the matrix, or the reverse). Retry the
        // exact/symbol tiers under each alternative symbol before falling back to the scan.
        let sym = core.rsplit('_').next().unwrap_or(core);
        for alias in alias_candidates(sym) {
            if let Some(&i) = self.exact.get(&alias).or_else(|| self.symbol.get(&alias)) {
                return Some(i);
            }
        }
        self.match_segments(&gl)
    }

    /// The fallback tier, a segment match over the lowered
    /// vocabulary: the first row that ends with `_{gl}`, starts with `{gl}_`
    /// or holds `_{gl}_`. In each case `gl` sits between `_` boundaries of
    /// the row, so the row holds `gl`'s first `_`-segment as a whole segment;
    /// only those rows are tested, in row order. The needles are built once
    /// and compared bytewise (the vocabulary and `gl` are both lowercase); an
    /// exact `row == gl` match belongs to the `exact` tier.
    fn match_segments(&self, gl: &str) -> Option<usize> {
        let first = gl.split('_').next().unwrap_or(gl);
        let rows = self.segment.get(first)?;
        let suffix = format!("_{gl}");
        let prefix = format!("{gl}_");
        let middle = format!("_{gl}_");
        rows.iter().copied().find(|&i| {
            let t = &self.lowered[i];
            t.ends_with(&suffix) || t.starts_with(&prefix) || t.contains(&middle)
        })
    }
}

/// Inverse-document-frequency marker weight `ln(C / df)`: a gene claimed by
/// all `C` types gets weight 0 (removed from scoring), a type-exclusive gene
/// the maximum `ln(C)`.
#[must_use]
pub fn idf_weight(n_types: usize, df: usize) -> f32 {
    (n_types as f32 / df.max(1) as f32).ln()
}

#[cfg(test)]
mod tests;
