//! Variant identifiers and genome builds.
//!
//! Association tables name a variant in one of a handful of spellings:
//! GTEx `chr1_12345_A_G_b38`, the eQTL Catalogue's `chr1_12345_A_G`, and
//! the colon form `1:12345:A:G` of OpenGWAS and many summary statistics.
//! [`parse_variant_id`] reads all of them; [`parse_locus`] accepts those or
//! any region / position [`crate::coordinates::parse_region`] reads, so a
//! column that mixes the two still resolves to coordinates.

use crate::coordinates::{chr_stripped, parse_region, PeakCoord};

/// A reference assembly of the human genome.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum GenomeBuild {
    GRCh37,
    GRCh38,
}

/// Length of chromosome 1 in each build: a VCF `##contig` line settles the
/// build even when it names no assembly.
const CHR1_LEN_GRCH37: u64 = 249_250_621;
const CHR1_LEN_GRCH38: u64 = 248_956_422;

impl GenomeBuild {
    /// Parse a build name: `GRCh38`, `hg38`, `b38`, `38` (and the 37/hg19
    /// equivalents), case-insensitive.
    pub fn parse(s: &str) -> Option<Self> {
        match s.trim().to_ascii_lowercase().as_str() {
            "grch38" | "hg38" | "b38" | "38" => Some(Self::GRCh38),
            "grch37" | "hg19" | "b37" | "37" => Some(Self::GRCh37),
            _ => None,
        }
    }

    pub fn as_str(&self) -> &'static str {
        match self {
            Self::GRCh37 => "GRCh37",
            Self::GRCh38 => "GRCh38",
        }
    }

    /// The build a piece of free text names — a file name, a VCF
    /// `##reference` / `##contig` line — or `None` when it names none or
    /// names both. Tokens are matched whole (`hg19`, `GRCh37`, `b37`,
    /// `hg38`, `GRCh38`, `b38`), so `b38` inside a longer word is ignored.
    pub fn detect(text: &str) -> Option<Self> {
        let lower = text.to_ascii_lowercase();
        let mut found: Option<Self> = None;
        for tok in lower.split(|c: char| !c.is_ascii_alphanumeric()) {
            let b = match tok {
                "grch38" | "hg38" | "b38" => Self::GRCh38,
                "grch37" | "hg19" | "b37" => Self::GRCh37,
                _ => continue,
            };
            match found {
                None => found = Some(b),
                Some(f) if f != b => return None,
                _ => {}
            }
        }
        found
    }

    /// The build implied by the length of chromosome 1, if it is one of the
    /// two.
    pub fn from_chr1_length(len: u64) -> Option<Self> {
        match len {
            CHR1_LEN_GRCH37 => Some(Self::GRCh37),
            CHR1_LEN_GRCH38 => Some(Self::GRCh38),
            _ => None,
        }
    }
}

impl std::fmt::Display for GenomeBuild {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

impl std::str::FromStr for GenomeBuild {
    type Err = anyhow::Error;
    fn from_str(s: &str) -> anyhow::Result<Self> {
        Self::parse(s)
            .ok_or_else(|| anyhow::anyhow!("unknown genome build `{s}` (GRCh37 or GRCh38)"))
    }
}

/// A variant: chromosome (no `chr` prefix), 1-based position, alleles, and
/// the build when the id carries one (`_b38`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VariantId {
    pub chr: Box<str>,
    pub pos: i64,
    pub ref_allele: Box<str>,
    pub alt_allele: Box<str>,
    pub build: Option<GenomeBuild>,
}

impl VariantId {
    /// The one-base locus `[pos, pos + 1)` of the variant.
    pub fn locus(&self) -> PeakCoord {
        PeakCoord {
            chr: self.chr.clone(),
            start: self.pos,
            end: self.pos + 1,
        }
    }
}

fn is_allele(s: &str) -> bool {
    !s.is_empty()
        && s.bytes().all(|b| {
            matches!(
                b.to_ascii_uppercase(),
                b'A' | b'C' | b'G' | b'T' | b'N' | b'*'
            )
        })
}

/// Parse `chr_pos_ref_alt[_b37|_b38]` or `chr:pos:ref:alt` (either with or
/// without the `chr` prefix). Alleles must be nucleotide strings; `None`
/// otherwise, so a name like `chr1_100_200` is not mistaken for one.
pub fn parse_variant_id(name: &str) -> Option<VariantId> {
    let name = name.trim();
    let sep = if name.contains(':') { ':' } else { '_' };
    let parts: Vec<&str> = name.split(sep).collect();
    let (chr, pos, r, a, build) = match parts.as_slice() {
        [c, p, r, a] => (*c, *p, *r, *a, None),
        [c, p, r, a, b] if sep == '_' => (*c, *p, *r, *a, Some(GenomeBuild::parse(b)?)),
        _ => return None,
    };
    if chr.is_empty() || !is_allele(r) || !is_allele(a) {
        return None;
    }
    let pos: i64 = pos.parse().ok()?;
    (pos > 0).then(|| VariantId {
        chr: chr_stripped(chr).into(),
        pos,
        ref_allele: r.into(),
        alt_allele: a.into(),
        build,
    })
}

/// Coordinates of any locus spelling: a variant id ([`parse_variant_id`]),
/// or a region / position ([`parse_region`]). The `chr` prefix is dropped.
pub fn parse_locus(name: &str) -> Option<PeakCoord> {
    parse_variant_id(name)
        .map(|v| v.locus())
        .or_else(|| parse_region(name))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variant_ids_parse_in_the_gtex_catalogue_and_colon_spellings() {
        let v = parse_variant_id("chr1_12345_A_G_b38").unwrap();
        assert_eq!((v.chr.as_ref(), v.pos), ("1", 12345));
        assert_eq!((v.ref_allele.as_ref(), v.alt_allele.as_ref()), ("A", "G"));
        assert_eq!(v.build, Some(GenomeBuild::GRCh38));
        let v = parse_variant_id("chrX_500_AT_A").unwrap();
        assert_eq!((v.chr.as_ref(), v.pos, v.build), ("X", 500, None));
        let v = parse_variant_id("1:12345:C:T").unwrap();
        assert_eq!((v.chr.as_ref(), v.pos), ("1", 12345));
        assert_eq!(parse_variant_id("chr2:7:G:GA").unwrap().chr.as_ref(), "2");
        assert!(parse_variant_id("chr1_100_200").is_none(), "an interval");
        assert!(parse_variant_id("chr1_100_A_G_b99").is_none());
        assert!(parse_variant_id("rs12345").is_none());
        assert!(parse_variant_id("1:100:A:G:b38").is_none());
    }

    #[test]
    fn loci_accept_variants_regions_and_positions() {
        assert_eq!(
            parse_locus("chr1_12345_A_G_b38").unwrap().to_string(),
            "1:12345-12346"
        );
        assert_eq!(
            parse_locus("1:12345:A:G").unwrap().to_string(),
            "1:12345-12346"
        );
        assert_eq!(
            parse_locus("chr1:100-200").unwrap().to_string(),
            "1:100-200"
        );
        assert!(
            parse_locus("chr1_100_200").is_none(),
            "regions are colon form"
        );
        assert_eq!(parse_locus("chr7:55").unwrap().to_string(), "7:55-56");
        assert!(parse_locus("GENE1").is_none());
    }

    #[test]
    fn builds_parse_and_detect_from_text() {
        assert_eq!(GenomeBuild::parse("hg19"), Some(GenomeBuild::GRCh37));
        assert_eq!(
            "GRCh38".parse::<GenomeBuild>().unwrap(),
            GenomeBuild::GRCh38
        );
        assert!("hg18".parse::<GenomeBuild>().is_err());
        assert_eq!(
            GenomeBuild::detect("ABC_hg19_K562.tsv.gz"),
            Some(GenomeBuild::GRCh37)
        );
        assert_eq!(
            GenomeBuild::detect("##contig=<ID=1,length=248956422,assembly=GRCh38>"),
            Some(GenomeBuild::GRCh38)
        );
        assert_eq!(
            GenomeBuild::detect("HG19/GRCh37"),
            Some(GenomeBuild::GRCh37)
        );
        assert_eq!(
            GenomeBuild::detect("liftover_hg19_to_hg38"),
            None,
            "both named"
        );
        assert_eq!(GenomeBuild::detect("rb38x.tsv"), None, "not a whole token");
        assert_eq!(
            GenomeBuild::from_chr1_length(249_250_621),
            Some(GenomeBuild::GRCh37)
        );
    }
}
