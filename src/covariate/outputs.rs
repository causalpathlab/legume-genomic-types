//! Tables a covariate fit writes beside its run, the same in every tool.

use super::{bin_of, CovariateModel};
use crate::batch_label;
use legume_numeric::matrix::parquet::{write_named_table, Column};
use nalgebra::DMatrix;
use std::collections::BTreeMap;

/// Each cell's scores (one column per functional, by its label) and
/// `batch_label`, rows named by `cells`, to `path` (parquet).
///
/// # Errors
/// A table that cannot be written.
pub fn write_cells(
    path: &std::path::Path,
    model: &CovariateModel,
    cells: &[Box<str>],
    scores: &DMatrix<f32>,
    labels: &[Box<str>],
) -> anyhow::Result<()> {
    let cols: Vec<Vec<f32>> = scores
        .column_iter()
        .map(|c| c.iter().copied().collect())
        .collect();
    let mut table: Vec<(Box<str>, Column)> = model
        .functionals
        .iter()
        .zip(&cols)
        .map(|(f, c)| (f.label.clone(), Column::F32(c)))
        .collect();
    table.push(("batch_label".into(), Column::Str(labels)));
    write_named_table(&path.to_string_lossy(), "cell", cells, &table)
}

/// Each gene's covariate effect, from the per-gene δ over the batch levels
/// (`log_delta`, genes × levels, the levels named): per bin, the mean over
/// batches of log δ at `batch|cvK`, centred over the bins, so a gene's row
/// says how much each bin was absorbed. Levels without a bin (plain labels)
/// are left out; nothing is written when no level has one. To `path`
/// (parquet); returns whether it was written.
///
/// # Errors
/// A table that cannot be written.
pub fn write_covariate_effect(
    path: &std::path::Path,
    genes: &[Box<str>],
    levels: &[Box<str>],
    log_delta: &DMatrix<f32>,
) -> anyhow::Result<bool> {
    let mut of_bin: BTreeMap<usize, (Box<str>, Vec<usize>)> = BTreeMap::new();
    for (l, level) in levels.iter().enumerate() {
        if let [_, bin] = batch_label::split(level).as_slice() {
            if let Some(k) = bin_of(bin) {
                of_bin
                    .entry(k)
                    .or_insert_with(|| ((*bin).into(), Vec::new()))
                    .1
                    .push(l);
            }
        }
    }
    if of_bin.is_empty() {
        return Ok(false);
    }
    let (bins, members): (Vec<Box<str>>, Vec<Vec<usize>>) = of_bin.into_values().unzip();
    let mut effect: Vec<Vec<f32>> = vec![vec![0.0; genes.len()]; bins.len()];
    for g in 0..genes.len() {
        let means: Vec<f64> = members
            .iter()
            .map(|ls| {
                ls.iter()
                    .map(|&l| f64::from(log_delta[(g, l)]))
                    .sum::<f64>()
                    / ls.len() as f64
            })
            .collect();
        let centre = means.iter().sum::<f64>() / means.len() as f64;
        for (k, m) in means.iter().enumerate() {
            effect[k][g] = (m - centre) as f32;
        }
    }
    let cols: Vec<(Box<str>, Column)> = bins
        .iter()
        .zip(&effect)
        .map(|(b, e)| (b.clone(), Column::F32(e)))
        .collect();
    write_named_table(&path.to_string_lossy(), "gene", genes, &cols)?;
    Ok(true)
}
