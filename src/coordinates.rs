use crate::gff::GeneId;
use crate::sam::Strand;
use rustc_hash::FxHashMap;

/// Compare chromosome names ignoring an optional "chr" prefix.
///
/// Handles mixed conventions (e.g., "chr1" vs "1", "chrX" vs "X"). The rest
/// of the name is compared as written, so "X" and "x" differ.
pub fn chr_eq(a: &str, b: &str) -> bool {
    chr_stripped(a) == chr_stripped(b)
}

/// Strip a "chr" prefix, in any case, from a chromosome name for use as a
/// lookup key. The rest keeps its case.
pub fn chr_stripped(s: &str) -> &str {
    match s.get(..3) {
        Some(p) if p.eq_ignore_ascii_case("chr") => &s[3..],
        _ => s,
    }
}

/// Parsed genomic coordinate for an ATAC peak.
#[derive(Debug, Clone)]
pub struct PeakCoord {
    pub chr: Box<str>,
    pub start: i64,
    pub end: i64,
}

impl PeakCoord {
    /// Canonical row key `{chr}:{start}-{end}`: "chr" prefix dropped, case
    /// kept (`chrX:0-100` becomes `X:0-100`). `chrX:0-100` and `X:0-100` get
    /// the same key, and the key parses back through [`parse_interval`] to
    /// the same interval. (It is its own key unless the chromosome, once
    /// stripped, starts with "chr" again, as in `chrChr1`.)
    pub fn locus_key(&self) -> Box<str> {
        format!("{}:{}-{}", chr_stripped(&self.chr), self.start, self.end).into_boxed_str()
    }
}

/// Gene TSS position parsed from GFF.
#[derive(Debug, Clone)]
pub struct GeneTss {
    pub chr: Box<str>,
    pub tss: i64,
}

/// Gene TSS position **with strand** parsed from GFF. Like [`GeneTss`]
/// but retains the strand so callers (e.g. strand-resolved genomic
/// pileups) can split forward (Watson) from backward (Crick) genes.
#[derive(Debug, Clone)]
pub struct GeneLoc {
    pub chr: Box<str>,
    /// Transcription start site (start on `+`, stop on `-`).
    pub tss: i64,
    pub strand: Strand,
}

/// Gene annotation simplified for eQTL and cis-regulatory analysis.
#[derive(Debug, Clone)]
pub struct Gene {
    pub gene_id: GeneId,
    pub gene_name: Option<Box<str>>,
    pub chromosome: Box<str>,
    pub tss: u64,
    pub strand: Strand,
}

/// Collection of gene annotations with cis window utilities.
#[derive(Debug, Clone)]
pub struct GeneAnnotations {
    pub genes: Vec<Gene>,
    pub cis_window: u64,
}

impl GeneAnnotations {
    /// Get cis-regulatory window for a gene (TSS ± cis_window).
    pub fn cis_region(&self, gene_idx: usize) -> (u64, u64) {
        let gene = &self.genes[gene_idx];
        let start = gene.tss.saturating_sub(self.cis_window);
        let end = gene.tss + self.cis_window;
        (start, end)
    }

    /// Get SNP indices within cis window for a gene.
    pub fn cis_snp_indices(
        &self,
        gene_idx: usize,
        snp_positions: &[u64],
        snp_chromosomes: &[Box<str>],
    ) -> Vec<usize> {
        let gene = &self.genes[gene_idx];
        let (start, end) = self.cis_region(gene_idx);

        snp_positions
            .iter()
            .enumerate()
            .filter(|(idx, &pos)| {
                chr_eq(&snp_chromosomes[*idx], &gene.chromosome) && pos >= start && pos <= end
            })
            .map(|(idx, _)| idx)
            .collect()
    }

    /// Filter genes to a specific genomic region.
    pub fn filter_to_region(&self, chromosome: &str, start: u64, end: u64) -> GeneAnnotations {
        let filtered_genes: Vec<Gene> = self
            .genes
            .iter()
            .filter(|g| chr_eq(g.chromosome.as_ref(), chromosome) && g.tss >= start && g.tss <= end)
            .cloned()
            .collect();

        GeneAnnotations {
            genes: filtered_genes,
            cis_window: self.cis_window,
        }
    }
}

/// The one coordinate grammar, colon form only: `chr:start-end`, and (when
/// `allow_position`) a single position `chr:pos` as the one-base interval
/// `[pos, pos + 1)`. The name holds exactly one `:`, so a pair id like
/// `chr1:1-2:3-4` is not one locus; the coordinates are plain digits (no
/// sign) and the name holds no whitespace. Underscore or dash spellings
/// (`chr1_100_200`, `chr1-100-200`) are not coordinates; see
/// [`import_interval`]. The chromosome comes back as written. `None` when
/// the name is not a coordinate.
fn parse_coordinate(name: &str, allow_position: bool) -> Option<PeakCoord> {
    let (chr, start, end) = coordinate_parts(name, allow_position)?;
    Some(PeakCoord {
        chr: chr.into(),
        start,
        end,
    })
}

/// [`parse_coordinate`] without allocating: the chromosome borrowed from
/// `name`.
fn coordinate_parts(name: &str, allow_position: bool) -> Option<(&str, i64, i64)> {
    let (chr, range) = name.split_once(':')?;
    match range.split_once('-') {
        Some((start, end)) => interval_parts(chr, start, end),
        None if allow_position => {
            let pos = plain_number(range)?;
            interval_parts_num(chr, pos, pos.checked_add(1)?)
        }
        None => None,
    }
}

/// A chromosome and two coordinate strings checked as one interval.
fn interval_parts<'a>(chr: &'a str, start: &str, end: &str) -> Option<(&'a str, i64, i64)> {
    interval_parts_num(chr, plain_number(start)?, plain_number(end)?)
}

fn interval_parts_num(chr: &str, start: i64, end: i64) -> Option<(&str, i64, i64)> {
    let chr_ok =
        !chr_stripped(chr).is_empty() && !chr.contains(|c: char| c == ':' || c.is_whitespace());
    (chr_ok && end > start).then_some((chr, start, end))
}

/// A coordinate: plain digits, no sign.
fn plain_number(s: &str) -> Option<i64> {
    s.parse::<i64>()
        .ok()
        .filter(|_| s.bytes().all(|b| b.is_ascii_digit()))
}

/// [`parse_interval`] without allocating: `(chr, start, end)`, the
/// chromosome borrowed from `name` as written.
pub fn split_interval(name: &str) -> Option<(&str, i64, i64)> {
    coordinate_parts(name, false)
}

/// True when `name` is a locus, `chr:start-end`; allocates nothing.
pub fn is_locus(name: &str) -> bool {
    split_interval(name).is_some()
}

/// The canonical key of a locus name (see [`PeakCoord::locus_key`]), or
/// `None` when `name` is not a locus. `chrX:0-100`, `CHRX:0-100` and
/// `X:0-100` all give `X:0-100`.
pub fn locus_key(name: &str) -> Option<Box<str>> {
    parse_interval(name).map(|l| l.locus_key())
}

/// Parse one interval name, `chr:start-end`, the chromosome kept verbatim.
/// Positions are not intervals; [`parse_region`] reads those, and
/// [`crate::variant::parse_locus`] also takes variant ids. `None` when the
/// name is not an interval.
pub fn parse_interval(name: &str) -> Option<PeakCoord> {
    parse_coordinate(name, false)
}

/// Read an interval as an outside tool spells it, for an importer only:
/// `chr:start-end`, `chr:start_end`, `chr-start-end` or `chr_start_end`
/// (Signac, some 10x and BED-derived exports). Without a `:` the numbers
/// are read from the right, so contig names with `_` or `-` stay whole.
/// The result's `Display` is the colon form every other function here
/// expects. Call it once, where the names enter, on rows known to be peaks;
/// a gene id like `GENE_1_2` would read as an interval too.
pub fn import_interval(name: &str) -> Option<PeakCoord> {
    let (chr, start, end) = match name.split_once(':') {
        Some((chr, range)) => {
            let (start, end) = range.split_once(['-', '_'])?;
            (chr, start, end)
        }
        None => {
            let mut parts = name.rsplitn(3, ['_', '-']);
            let (end, start) = (parts.next()?, parts.next()?);
            (parts.next()?, start, end)
        }
    };
    let (chr, start, end) = interval_parts(chr, start, end)?;
    Some(PeakCoord {
        chr: chr.into(),
        start,
        end,
    })
}

/// [`parse_interval`] over a list of peak names.
pub fn parse_peak_coordinates(peak_names: &[Box<str>]) -> Vec<Option<PeakCoord>> {
    peak_names.iter().map(|name| parse_interval(name)).collect()
}

/// Parse one genomic region name: an interval `chr:start-end`, or a single
/// position `chr:pos` (a SNP), which becomes the one-base interval
/// `[pos, pos + 1)`. The `chr` prefix is dropped so `chr1` and `1` name the
/// same chromosome. `None` when the name is not a coordinate.
pub fn parse_region(name: &str) -> Option<PeakCoord> {
    let mut r = parse_coordinate(name, true)?;
    r.chr = chr_stripped(&r.chr).into();
    Some(r)
}

/// Tile a region onto fixed windows `[i·w, (i+1)·w)`: every window the
/// region overlaps, ascending, on the region's chromosome. A window of 0 is
/// the region itself.
pub fn tile_windows(region: &PeakCoord, window: i64) -> Vec<PeakCoord> {
    if window <= 0 {
        return vec![region.clone()];
    }
    let first = region.start.div_euclid(window);
    let last = (region.end - 1).max(region.start).div_euclid(window);
    (first..=last)
        .map(|i| PeakCoord {
            chr: region.chr.clone(),
            start: i * window,
            end: (i + 1) * window,
        })
        .collect()
}

impl std::fmt::Display for PeakCoord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}-{}", self.chr, self.start, self.end)
    }
}

/// Find peaks within a cis window of a gene's TSS.
pub fn find_cis_peaks(
    gene_tss: &GeneTss,
    peak_coords: &[Option<PeakCoord>],
    window: i64,
) -> Vec<usize> {
    peak_coords
        .iter()
        .enumerate()
        .filter_map(|(idx, coord)| {
            let coord = coord.as_ref()?;
            if !chr_eq(coord.chr.as_ref(), gene_tss.chr.as_ref()) {
                return None;
            }
            let mid = (coord.start + coord.end) / 2;
            if (mid - gene_tss.tss).abs() <= window {
                Some(idx)
            } else {
                None
            }
        })
        .collect()
}

/// Load gene TSS positions from a GFF/GTF file, aligned to `gene_names`
/// via [`align_gene_loci`] (Ensembl id / HGNC symbol /
/// `ENSG…_SYMBOL`). Strand-discarding projection of [`load_gene_loci`].
pub fn load_gene_tss(
    gff_file: &str,
    gene_names: &[Box<str>],
) -> anyhow::Result<Vec<Option<GeneTss>>> {
    Ok(load_gene_loci(gff_file, gene_names)?
        .into_iter()
        .map(|loc| {
            loc.map(|l| GeneTss {
                chr: l.chr,
                tss: l.tss,
            })
        })
        .collect())
}

/// Load a `gene-symbol / Ensembl-id → (chr, TSS, strand)` map from a GFF/GTF.
///
/// Only `gene` features are kept. Each gene is registered under its HGNC
/// symbol (when present) **and** its version-stripped Ensembl id, so
/// [`align_gene_loci`] / [`GeneIndexResolver`] can match either form (and
/// compound `ENSG…_SYMBOL` row names).
pub fn load_gene_loci_map(gff_file: &str) -> anyhow::Result<FxHashMap<Box<str>, GeneLoc>> {
    use crate::gff::{parse_ensembl_id, read_gff_record_vec, FeatureType, GeneId, GeneSymbol};
    use log::info;

    let records = read_gff_record_vec(gff_file)?;

    let mut loc_map: FxHashMap<Box<str>, GeneLoc> = FxHashMap::default();
    for rec in &records {
        if rec.feature_type != FeatureType::Gene {
            continue;
        }
        let tss = match rec.strand {
            Strand::Forward => rec.start,
            Strand::Backward => rec.stop,
        };
        let loc = GeneLoc {
            chr: rec.seqname.clone(),
            tss,
            strand: rec.strand,
        };
        if let GeneSymbol::Symbol(s) = &rec.gene_name {
            loc_map.entry(s.clone()).or_insert_with(|| loc.clone());
        }
        if let GeneId::Ensembl(id) = &rec.gene_id {
            let key: Box<str> = parse_ensembl_id(id).unwrap_or(id.as_ref()).into();
            loc_map.entry(key).or_insert(loc);
        } else if matches!(rec.gene_name, GeneSymbol::Missing) {
            let id: Box<str> = rec.gene_id.clone().into();
            loc_map.entry(id).or_insert(loc);
        }
    }

    info!("Loaded {} gene-locus keys from GFF", loc_map.len());
    Ok(loc_map)
}

/// Align a `key → GeneLoc` map onto `gene_names` with
/// [`legume_numeric::matrix::membership::GeneIndexResolver`]: exact name,
/// HGNC symbol, or Ensembl id (including `ENSG…_SYMBOL` compounds).
pub fn align_gene_loci(
    gene_names: &[Box<str>],
    loci_by_key: &FxHashMap<Box<str>, GeneLoc>,
) -> Vec<Option<GeneLoc>> {
    align_by_gene_key(gene_names, loci_by_key)
}

/// Align a `key → GeneTss` map the same way as [`align_gene_loci`].
pub fn align_gene_tss(
    gene_names: &[Box<str>],
    tss_by_key: &FxHashMap<Box<str>, GeneTss>,
) -> Vec<Option<GeneTss>> {
    align_by_gene_key(gene_names, tss_by_key)
}

/// `values_by_key` onto `gene_names` by [`GeneIndexResolver`]; the first key
/// that resolves to a gene wins.
fn align_by_gene_key<T: Clone>(
    gene_names: &[Box<str>],
    values_by_key: &FxHashMap<Box<str>, T>,
) -> Vec<Option<T>> {
    use legume_numeric::matrix::membership::GeneIndexResolver;

    let resolver = GeneIndexResolver::build(gene_names, Some('_'), false);
    let mut out = vec![None; gene_names.len()];
    for (key, value) in values_by_key {
        if let Some(i) = resolver.resolve(key) {
            if out[i].is_none() {
                out[i] = Some(value.clone());
            }
        }
    }
    out
}

/// Load gene TSS positions **and strand** from a GFF/GTF file, aligned
/// to `gene_names` via [`align_gene_loci`] (`None` where unmatched).
pub fn load_gene_loci(
    gff_file: &str,
    gene_names: &[Box<str>],
) -> anyhow::Result<Vec<Option<GeneLoc>>> {
    use log::info;

    let loc_map = load_gene_loci_map(gff_file)?;
    let result = align_gene_loci(gene_names, &loc_map);

    let matched = result.iter().filter(|x| x.is_some()).count();
    info!("Matched {}/{} genes to GFF loci", matched, gene_names.len());

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regions_parse_intervals_and_positions_without_the_chr_prefix() {
        let r = parse_region("chr1:1000-2000").unwrap();
        assert_eq!((r.chr.as_ref(), r.start, r.end), ("1", 1000, 2000));
        let r = parse_region("chrX:5000").unwrap();
        assert_eq!((r.chr.as_ref(), r.start, r.end), ("X", 5000, 5001));
        assert!(parse_region("GENE1").is_none());
        assert!(parse_region("chr1:2000-1000").is_none(), "empty interval");
        assert!(parse_region(":1-2").is_none());
        assert_eq!(parse_region("chr2:10-20").unwrap().to_string(), "2:10-20");
    }

    #[test]
    fn only_the_colon_form_is_a_coordinate() {
        for name in [
            "1_1000_2000",
            "chr1_1000_2000",
            "chr1-1000-2000",
            "X_5000",
            "chr1:1000_2000",
            "chr1:+1000-2000",
            "chr1:1000-2000 ",
            " chr1:1000-2000",
            "chr1:1000-2000:+",
        ] {
            assert!(parse_region(name).is_none(), "{name}");
        }
    }

    #[test]
    fn loci_keep_chromosome_case_and_contig_names() {
        let l = parse_interval("chrX:0-100").unwrap();
        assert_eq!((l.chr.as_ref(), l.start, l.end), ("chrX", 0, 100));
        assert_eq!(l.locus_key().as_ref(), "X:0-100");
        for name in ["chrX:0-100", "X:0-100", "CHRX:0-100"] {
            assert_eq!(
                parse_interval(name).unwrap().locus_key().as_ref(),
                "X:0-100"
            );
        }
        let l = parse_interval("chrUn_CTG1v1:0-100").unwrap();
        assert_eq!(l.chr.as_ref(), "chrUn_CTG1v1");
        let l = parse_interval("chrUn_CTG1v1:0-100").unwrap();
        let key = l.locus_key();
        assert_eq!(
            parse_interval(&key).unwrap().locus_key(),
            key,
            "key round-trips"
        );
        // Exactly one `:`: a pair id or a name carrying `:` is not one locus.
        assert!(parse_interval("chr1:100-200:300-400").is_none());
        assert!(parse_interval("chr1:1-2:chr1:3-4").is_none());
        assert!(parse_interval("GENE1*01:01:01:5-10").is_none());
        assert!(
            parse_region("chr1:9223372036854775807").is_none(),
            "no overflow"
        );
        assert!(
            parse_interval("chr1:5000").is_none(),
            "a position is not an interval"
        );
        assert!(parse_interval("ENSG000_GENE1").is_none());
        assert!(parse_interval("GENE_1_2").is_none());
        assert!(parse_interval("GENE1*01:01:01").is_none());
        assert!(parse_interval("chr:1-2").is_none(), "no chromosome left");
        assert_eq!(locus_key("CHRX:0-100").as_deref(), Some("X:0-100"));
        assert!(is_locus("X:0-100") && !is_locus("X_0_100"));
        assert_eq!(split_interval("chrX:0-100"), Some(("chrX", 0, 100)));
        for name in [
            "chr1-100-200",
            "chr1_100_200",
            "chr1:100-200",
            "chr1:100_200",
        ] {
            assert_eq!(import_interval(name).unwrap().to_string(), "chr1:100-200");
        }
        let l = import_interval("chrUn_CTG1v1-5-10").unwrap();
        assert_eq!(l.to_string(), "chrUn_CTG1v1:5-10");
        assert!(import_interval("GENE1").is_none());
        assert!(chr_eq("chrX", "X") && chr_eq("cHrX", "X") && !chr_eq("x", "X"));
    }

    #[test]
    fn tiling_covers_every_overlapped_window_and_only_those() {
        let r = parse_region("chr1:4999-10001").unwrap();
        let w: Vec<String> = tile_windows(&r, 5000)
            .iter()
            .map(ToString::to_string)
            .collect();
        assert_eq!(w, vec!["1:0-5000", "1:5000-10000", "1:10000-15000"]);
        let snp = parse_region("chr1:5000").unwrap();
        let w: Vec<String> = tile_windows(&snp, 5000)
            .iter()
            .map(ToString::to_string)
            .collect();
        assert_eq!(
            w,
            vec!["1:5000-10000"],
            "a position at a boundary lands in one window"
        );
        let exact = parse_region("chr1:5000-10000").unwrap();
        assert_eq!(
            tile_windows(&exact, 5000).len(),
            1,
            "an exact window is one window"
        );
        assert_eq!(
            tile_windows(&exact, 0).len(),
            1,
            "window 0 keeps the region"
        );
        assert_eq!(tile_windows(&exact, 0)[0].to_string(), "1:5000-10000");
    }

    #[test]
    fn load_gene_loci_keeps_strand_and_tss() {
        // Two genes: AAA on +, BBB on -. TSS = start on +, stop on -.
        let gtf = "\
chr1\tHAVANA\tgene\t100\t200\t.\t+\t.\tgene_id \"ENSG001\"; gene_name \"AAA\"; gene_type \"protein_coding\"
chr1\tHAVANA\tgene\t400\t600\t.\t-\t.\tgene_id \"ENSG002\"; gene_name \"BBB\"; gene_type \"protein_coding\"
chr1\tHAVANA\texon\t100\t150\t.\t+\t.\tgene_id \"ENSG001\"; gene_name \"AAA\"
";
        let path = std::env::temp_dir().join(format!("genloci_{}.gtf", std::process::id()));
        std::fs::write(&path, gtf).unwrap();
        let names: Vec<Box<str>> = vec!["AAA".into(), "BBB".into(), "MISSING".into()];
        let loci = load_gene_loci(path.to_str().unwrap(), &names).unwrap();
        std::fs::remove_file(&path).ok();

        let aaa = loci[0].as_ref().expect("AAA present");
        assert!(matches!(aaa.strand, Strand::Forward));
        assert_eq!(aaa.tss, 100);
        assert_eq!(chr_stripped(&aaa.chr), "1");

        let bbb = loci[1].as_ref().expect("BBB present");
        assert!(matches!(bbb.strand, Strand::Backward));
        assert_eq!(bbb.tss, 600);

        assert!(loci[2].is_none());
    }

    #[test]
    fn load_gene_loci_matches_ensg_symbol_compounds() {
        let gtf = "\
chr1\tHAVANA\tgene\t100\t200\t.\t+\t.\tgene_id \"ENSG00000186092.7\"; gene_name \"OR4F5\"; gene_type \"protein_coding\"
chr1\tHAVANA\tgene\t400\t600\t.\t-\t.\tgene_id \"ENSG00000237613.2\"; gene_name \"FAM138A\"; gene_type \"protein_coding\"
";
        let path = std::env::temp_dir().join(format!("genloci_cmp_{}.gtf", std::process::id()));
        std::fs::write(&path, gtf).unwrap();
        let names: Vec<Box<str>> = vec![
            "ENSG00000186092_OR4F5".into(),
            "FAM138A".into(),
            "ENSG00000237613".into(),
        ];
        let loci = load_gene_loci(path.to_str().unwrap(), &names).unwrap();
        std::fs::remove_file(&path).ok();
        assert_eq!(loci[0].as_ref().unwrap().tss, 100);
        assert_eq!(loci[1].as_ref().unwrap().tss, 600);
        assert_eq!(loci[2].as_ref().unwrap().tss, 600);
    }
}
