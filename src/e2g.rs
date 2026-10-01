//! The E2G enhancer–gene release (`gs://e2g`): normalized parquet tables
//! joined back into one enhancer → gene link per prediction.
//!
//! The release is a directory of three parts: `cell_types.parquet` (`id`,
//! `name`), `enhancers.parquet` (`id`, `chromosome`, `start`, `end`; BED,
//! 0-based) and `enhancer_gene_predictions/*.parquet`, one file per
//! chromosome (`enhancer_id`, target gene, `score`, `model`,
//! `cell_type_id`). [`E2gRelease::open`] loads the two small tables;
//! [`E2gRelease::for_each_link`] streams the predictions and hands every one
//! to a callback with its enhancer's coordinates and its cell type's name
//! resolved, borrowing the row rather than allocating per link.

use legume_numeric::matrix::table::TableReader;
use rustc_hash::FxHashMap;
use std::path::{Path, PathBuf};

/// `cell_types.parquet`: the file that marks a directory as a release.
pub const CELL_TYPES: &str = "cell_types.parquet";
const ENHANCERS: &str = "enhancers.parquet";
const PREDICTIONS: &str = "enhancer_gene_predictions";

/// One prediction, borrowed from the row being read.
#[derive(Debug)]
pub struct E2gLink<'a> {
    /// Chromosome as the prediction file writes it (`17`, `X`).
    pub chr: &'a str,
    /// The enhancer as a 1-based half-open `[start, end)`; `None` when its
    /// id is not in `enhancers.parquet`.
    pub span: Option<(i64, i64)>,
    /// Target gene symbol, or its Ensembl id when the symbol is empty.
    pub gene: &'a str,
    /// `None` when the score is not a number.
    pub score: Option<f64>,
    /// `ENCODE-rE2G`, `scE2G`, ….
    pub model: &'a str,
    /// The sample's cell type name (`cell_types.parquet` `name`, treatment
    /// included); `None` when the sample id is not listed.
    pub cell_type: Option<&'a str>,
}

/// A release opened for streaming.
pub struct E2gRelease {
    /// Sample id → cell type name.
    cell_types: FxHashMap<Box<str>, Box<str>>,
    /// Enhancer id → 1-based `[start, end)`; `(0, 0)` for ids not present.
    spans: Vec<(i64, i64)>,
    files: Vec<PathBuf>,
}

impl E2gRelease {
    /// `true` when `dir` holds a release (it has `cell_types.parquet`).
    pub fn is_release(dir: &Path) -> bool {
        dir.join(CELL_TYPES).exists()
    }

    /// Load the cell types and the enhancer coordinates, and list the
    /// prediction files (sorted).
    pub fn open(dir: &Path) -> anyhow::Result<Self> {
        let path = |f: &str| dir.join(f).to_string_lossy().into_owned();

        let t = TableReader::open(&path(CELL_TYPES))?;
        let cols = t.select(&["id", "name"])?;
        let mut cell_types = FxHashMap::default();
        for row in t.rows(&cols)? {
            let mut row = row?;
            let name = std::mem::take(&mut row[1]);
            cell_types.insert(std::mem::take(&mut row[0]), name);
        }

        let t = TableReader::open(&path(ENHANCERS))?;
        let cols = t.select(&["id", "start", "end"])?;
        let mut spans: Vec<(i64, i64)> = Vec::new();
        for row in t.rows(&cols)? {
            let row = row?;
            let (Ok(id), Ok(s), Ok(e)) = (
                row[0].parse::<usize>(),
                row[1].parse::<i64>(),
                row[2].parse::<i64>(),
            ) else {
                continue;
            };
            if spans.len() <= id {
                spans.resize(id + 1, (0, 0));
            }
            // BED is 0-based half-open; shift to 1-based.
            spans[id] = (s + 1, e + 1);
        }

        let mut files: Vec<PathBuf> = std::fs::read_dir(dir.join(PREDICTIONS))
            .map_err(|e| anyhow::anyhow!("{}: {e}", dir.join(PREDICTIONS).display()))?
            .filter_map(|e| e.ok().map(|e| e.path()))
            .filter(|p| p.extension().is_some_and(|x| x == "parquet"))
            .collect();
        files.sort();
        anyhow::ensure!(
            !files.is_empty(),
            "{}: no {PREDICTIONS}/*.parquet",
            dir.display()
        );
        Ok(Self {
            cell_types,
            spans,
            files,
        })
    }

    /// Number of samples in `cell_types.parquet`.
    pub fn n_samples(&self) -> usize {
        self.cell_types.len()
    }

    /// Stream every prediction, file by file, into `f`.
    pub fn for_each_link(&self, mut f: impl FnMut(E2gLink<'_>)) -> anyhow::Result<()> {
        for file in &self.files {
            let t = TableReader::open(&file.to_string_lossy())?;
            let cols = t.select(&[
                "chromosome",
                "enhancer_id",
                "target_gene_name",
                "target_gene_id",
                "score",
                "model",
                "cell_type_id",
            ])?;
            for row in t.rows(&cols)? {
                let row = row?;
                let span = row[1]
                    .parse::<usize>()
                    .ok()
                    .and_then(|id| self.spans.get(id).copied())
                    .filter(|&(s, e)| e > s);
                let gene = match row[2].trim() {
                    "" => row[3].trim(),
                    g => g,
                };
                f(E2gLink {
                    chr: row[0].trim(),
                    span,
                    gene,
                    score: row[4].trim().parse::<f64>().ok(),
                    model: row[5].trim(),
                    cell_type: self.cell_types.get(row[6].trim()).map(AsRef::as_ref),
                });
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use legume_numeric::matrix::parquet::{write_table, Column};

    fn s(v: &[&str]) -> Vec<Box<str>> {
        v.iter().map(|x| Box::from(*x)).collect()
    }

    #[test]
    fn predictions_come_back_joined_to_coordinates_and_cell_type_names() {
        let dir = tempfile::tempdir().unwrap();
        let root = dir.path();
        std::fs::create_dir(root.join(PREDICTIONS)).unwrap();
        let at = |f: &str| root.join(f).to_string_lossy().into_owned();
        write_table(
            &at(CELL_TYPES),
            &[
                ("id".into(), Column::Str(&s(&["HepG2_A", "HepG2_B"]))),
                ("name".into(), Column::Str(&s(&["HepG2", "HepG2"]))),
            ],
        )
        .unwrap();
        write_table(
            &at(ENHANCERS),
            &[
                ("id".into(), Column::I32(&[3])),
                ("chromosome".into(), Column::Str(&s(&["17"]))),
                ("start".into(), Column::I32(&[99])),
                ("end".into(), Column::I32(&[200])),
            ],
        )
        .unwrap();
        write_table(
            &at("enhancer_gene_predictions/chr17.parquet"),
            &[
                ("chromosome".into(), Column::Str(&s(&["17", "17"]))),
                ("enhancer_id".into(), Column::I32(&[3, 9])),
                ("target_gene_name".into(), Column::Str(&s(&["", "GENE1"]))),
                (
                    "target_gene_id".into(),
                    Column::Str(&s(&["ENSG1", "ENSG2"])),
                ),
                ("score".into(), Column::F32(&[0.5, 0.25])),
                ("model".into(), Column::Str(&s(&["ENCODE-rE2G", "scE2G"]))),
                ("cell_type_id".into(), Column::Str(&s(&["HepG2_B", "gone"]))),
            ],
        )
        .unwrap();
        assert!(E2gRelease::is_release(root));
        let r = E2gRelease::open(root).unwrap();
        assert_eq!(r.n_samples(), 2);
        let mut got = Vec::new();
        r.for_each_link(|l| {
            got.push((
                l.chr.to_string(),
                l.span,
                l.gene.to_string(),
                l.score,
                l.model.to_string(),
                l.cell_type.map(str::to_string),
            ))
        })
        .unwrap();
        assert_eq!(
            got,
            vec![
                (
                    "17".into(),
                    Some((100, 201)),
                    "ENSG1".into(),
                    Some(0.5),
                    "ENCODE-rE2G".into(),
                    Some("HepG2".into())
                ),
                (
                    "17".into(),
                    None,
                    "GENE1".into(),
                    Some(0.25),
                    "scE2G".into(),
                    None
                ),
            ]
        );
    }
}
