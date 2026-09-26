//! A streaming VCF reader for summary statistics and position lists.
//!
//! Two shapes are in scope. A *GWAS-VCF* (the MRC IEU / OpenGWAS format)
//! carries one sample column per trait with `FORMAT` keys `ES:SE:LP[:…]`:
//! effect size, its standard error, and `-log10 p`. A *position-only* VCF
//! (a list of fine-mapped or lead variants) has no samples at all. Both are
//! read one record at a time, gzip or not, so a 10-million-variant file is
//! never held in memory. Genotype VCFs are not the target: sample values
//! are kept as text, with no genotype parsing.

use crate::coordinates::chr_stripped;
use crate::variant::GenomeBuild;
use legume_numeric::matrix::common_io::open_buf_reader;
use std::io::BufRead;

/// The `##` meta lines and the column header of a VCF.
#[derive(Clone, Debug, Default)]
pub struct VcfHeader {
    /// Every `##key=value` line, `##` stripped, in file order.
    pub meta: Vec<Box<str>>,
    /// Sample (trait) column names after `FORMAT`.
    pub samples: Vec<Box<str>>,
    /// The build the meta lines name (`##reference`, `##contig …
    /// assembly=`, or the length of contig 1), if they agree on one.
    pub build: Option<GenomeBuild>,
}

impl VcfHeader {
    /// Values of the `##key=` lines, `key` matched exactly.
    pub fn meta_values<'a>(&'a self, key: &'a str) -> impl Iterator<Item = &'a str> + 'a {
        self.meta
            .iter()
            .filter_map(move |m| m.strip_prefix(key).and_then(|rest| rest.strip_prefix('=')))
    }

    /// One `key=value` field out of a structured meta value
    /// `<ID=x,key=value,…>` (quotes stripped).
    pub fn structured_field<'a>(value: &'a str, key: &str) -> Option<&'a str> {
        let inner = value.trim().strip_prefix('<')?.strip_suffix('>')?;
        inner.split(',').find_map(|kv| {
            let (k, v) = kv.split_once('=')?;
            (k.trim() == key).then(|| v.trim().trim_matches('"'))
        })
    }

    fn detect_build(meta: &[Box<str>]) -> Option<GenomeBuild> {
        let mut found: Option<GenomeBuild> = None;
        let mut vote = |b: Option<GenomeBuild>| -> bool {
            match (found, b) {
                (_, None) => true,
                (None, Some(b)) => {
                    found = Some(b);
                    true
                }
                (Some(f), Some(b)) => f == b,
            }
        };
        for m in meta {
            let ok = if let Some(v) = m.strip_prefix("reference=") {
                vote(GenomeBuild::detect(v))
            } else if let Some(v) = m.strip_prefix("contig=") {
                let by_name = Self::structured_field(v, "assembly").and_then(GenomeBuild::detect);
                let by_len = match Self::structured_field(v, "ID").map(chr_stripped) {
                    Some("1") => Self::structured_field(v, "length")
                        .and_then(|l| l.parse().ok())
                        .and_then(GenomeBuild::from_chr1_length),
                    _ => None,
                };
                vote(by_name) && vote(by_len)
            } else {
                true
            };
            if !ok {
                return None;
            }
        }
        found
    }
}

/// One VCF data line.
#[derive(Clone, Debug)]
pub struct VcfRecord {
    /// Chromosome, `chr` prefix dropped.
    pub chr: Box<str>,
    /// 1-based position.
    pub pos: i64,
    pub id: Box<str>,
    pub ref_allele: Box<str>,
    pub alt_allele: Box<str>,
    pub info: Box<str>,
    /// `FORMAT` keys (empty for a position-only VCF).
    pub format: Vec<Box<str>>,
    /// Per sample, the raw `:`-joined column; split on demand by
    /// [`VcfRecord::sample_value`], so a reader that wants two keys does not
    /// pay for all of them.
    pub samples: Vec<Box<str>>,
}

/// GWAS-VCF statistics of one record for one sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GwasStat {
    pub es: Option<f64>,
    pub se: Option<f64>,
    /// `-log10 p`.
    pub lp: Option<f64>,
}

impl VcfRecord {
    /// The value of `FORMAT` key `key` for sample `s`.
    pub fn sample_value(&self, s: usize, key: &str) -> Option<&str> {
        let k = self.format.iter().position(|f| f.as_ref() == key)?;
        let v = self.samples.get(s)?.split(':').nth(k)?;
        (v != "." && !v.is_empty()).then_some(v)
    }

    /// `ES`, `SE` and `LP` of sample `s`; `None` fields where absent or `.`.
    pub fn gwas(&self, s: usize) -> GwasStat {
        let num = |k| self.sample_value(s, k).and_then(|v| v.parse::<f64>().ok());
        GwasStat {
            es: num("ES"),
            se: num("SE"),
            lp: num("LP"),
        }
    }
}

/// A VCF opened for streaming: the header is read by [`VcfReader::open`],
/// the records by iterating.
pub struct VcfReader {
    path: Box<str>,
    header: VcfHeader,
    reader: Box<dyn BufRead>,
    line: String,
    line_no: usize,
}

impl VcfReader {
    pub fn open(path: &str) -> anyhow::Result<Self> {
        let mut reader =
            open_buf_reader(path).map_err(|e| anyhow::anyhow!("opening {path}: {e}"))?;
        let mut meta: Vec<Box<str>> = Vec::new();
        let mut samples = Vec::new();
        let mut line = String::new();
        let mut line_no = 0usize;
        let mut saw_columns = false;
        loop {
            line.clear();
            if reader.read_line(&mut line)? == 0 {
                break;
            }
            line_no += 1;
            let l = line.trim_end_matches(['\n', '\r']);
            if let Some(m) = l.strip_prefix("##") {
                meta.push(m.into());
            } else if let Some(cols) = l.strip_prefix('#') {
                let cols: Vec<&str> = cols.split('\t').collect();
                anyhow::ensure!(
                    cols.first().map(|c| c.trim()) == Some("CHROM"),
                    "{path}: line {line_no}: expected the `#CHROM` column header"
                );
                samples = cols.iter().skip(9).map(|s| Box::from(s.trim())).collect();
                saw_columns = true;
                break;
            } else {
                anyhow::bail!("{path}: line {line_no}: data before the `#CHROM` header");
            }
        }
        anyhow::ensure!(saw_columns, "{path}: no `#CHROM` header; not a VCF");
        anyhow::ensure!(
            meta.first()
                .is_some_and(|m| m.starts_with("fileformat=VCF")),
            "{path}: the first line is not `##fileformat=VCF…`"
        );
        let build = VcfHeader::detect_build(&meta);
        Ok(Self {
            path: path.into(),
            header: VcfHeader {
                meta,
                samples,
                build,
            },
            reader,
            line,
            line_no,
        })
    }

    pub fn header(&self) -> &VcfHeader {
        &self.header
    }

    /// `true` when the `FORMAT` definitions name `ES`, `SE` and `LP`: a
    /// GWAS-VCF rather than a genotype or position VCF.
    pub fn is_gwas_vcf(&self) -> bool {
        let ids: Vec<&str> = self
            .header
            .meta_values("FORMAT")
            .filter_map(|v| VcfHeader::structured_field(v, "ID"))
            .collect();
        ["ES", "SE", "LP"].iter().all(|k| ids.contains(k))
    }

    fn parse_line(&self, l: &str) -> anyhow::Result<VcfRecord> {
        let f: Vec<&str> = l.split('\t').collect();
        anyhow::ensure!(
            f.len() >= 5,
            "{}: line {}: {} columns, a VCF record needs at least 5",
            self.path,
            self.line_no,
            f.len()
        );
        let pos: i64 = f[1].trim().parse().map_err(|e| {
            anyhow::anyhow!("{}: line {}: POS `{}`: {e}", self.path, self.line_no, f[1])
        })?;
        let format: Vec<Box<str>> = match f.get(8) {
            Some(fmt) if !fmt.is_empty() && *fmt != "." => fmt.split(':').map(Box::from).collect(),
            _ => Vec::new(),
        };
        let samples = f.iter().skip(9).map(|s| Box::from(*s)).collect();
        Ok(VcfRecord {
            chr: chr_stripped(f[0].trim()).into(),
            pos,
            id: f[2].into(),
            ref_allele: f[3].into(),
            alt_allele: f[4].into(),
            info: f.get(7).copied().unwrap_or(".").into(),
            format,
            samples,
        })
    }
}

impl Iterator for VcfReader {
    type Item = anyhow::Result<VcfRecord>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.line.clear();
            match self.reader.read_line(&mut self.line) {
                Ok(0) => return None,
                Ok(_) => {}
                Err(e) => return Some(Err(e.into())),
            }
            self.line_no += 1;
            let l = self.line.trim_end_matches(['\n', '\r']);
            if l.is_empty() || l.starts_with('#') {
                continue;
            }
            // Parse out of a taken buffer: no copy of the line.
            let line = std::mem::take(&mut self.line);
            let rec = self.parse_line(line.trim_end_matches(['\n', '\r']));
            self.line = line;
            return Some(rec);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GWAS_VCF: &str = "##fileformat=VCFv4.2
##FORMAT=<ID=ES,Number=A,Type=Float,Description=\"Effect size estimate relative to the alternative allele\">
##FORMAT=<ID=SE,Number=A,Type=Float,Description=\"Standard error of effect size estimate\">
##FORMAT=<ID=LP,Number=A,Type=Float,Description=\"-log10 p-value for effect estimate\">
##FORMAT=<ID=ID,Number=1,Type=String,Description=\"Study variant identifier\">
##contig=<ID=1,length=249250621,assembly=HG19/GRCh37>
##SAMPLE=<ID=ieu-a-2,TotalVariants=2555511,StudyType=Continuous>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tieu-a-2
chr1\t100\trs1\tA\tG\t.\tPASS\t.\tES:SE:LP:ID\t0.05:0.01:8.5:rs1
1\t200\trs2\tC\tT\t.\tPASS\t.\tES:SE:LP:ID\t-0.01:0.02:.:rs2
";

    #[test]
    fn a_gwas_vcf_header_names_its_trait_its_build_and_its_statistics() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("ieu-a-2.vcf");
        std::fs::write(&p, GWAS_VCF).unwrap();
        let r = VcfReader::open(p.to_str().unwrap()).unwrap();
        assert_eq!(r.header().samples, vec![Box::<str>::from("ieu-a-2")]);
        assert_eq!(r.header().build, Some(GenomeBuild::GRCh37));
        assert!(r.is_gwas_vcf());
        let sample = r.header().meta_values("SAMPLE").next().unwrap();
        assert_eq!(
            VcfHeader::structured_field(sample, "StudyType"),
            Some("Continuous")
        );
        let recs: Vec<VcfRecord> = r.map(Result::unwrap).collect();
        assert_eq!(recs.len(), 2);
        assert_eq!((recs[0].chr.as_ref(), recs[0].pos), ("1", 100));
        assert_eq!(
            recs[0].gwas(0),
            GwasStat {
                es: Some(0.05),
                se: Some(0.01),
                lp: Some(8.5)
            }
        );
        assert_eq!(recs[1].gwas(0).lp, None, "`.` is missing");
        assert_eq!(recs[1].sample_value(0, "ID"), Some("rs2"));
    }

    #[test]
    fn a_position_only_vcf_has_no_samples_and_a_build_from_contig_length() {
        let body = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=248956422>\n\
                    #CHROM\tPOS\tID\tREF\tALT\n\
                    chr1\t5\t.\tA\tC\n";
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("pos.vcf");
        std::fs::write(&p, body).unwrap();
        let r = VcfReader::open(p.to_str().unwrap()).unwrap();
        assert!(r.header().samples.is_empty());
        assert!(!r.is_gwas_vcf());
        assert_eq!(r.header().build, Some(GenomeBuild::GRCh38));
        let recs: Vec<VcfRecord> = r.map(Result::unwrap).collect();
        assert_eq!(recs[0].pos, 5);
        assert!(recs[0].format.is_empty());
    }

    #[test]
    fn a_file_that_is_not_a_vcf_is_refused() {
        let dir = tempfile::tempdir().unwrap();
        let p = dir.path().join("x.vcf");
        std::fs::write(&p, "a\tb\n1\t2\n").unwrap();
        assert!(VcfReader::open(p.to_str().unwrap()).is_err());
    }
}
