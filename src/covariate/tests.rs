//! Covariates end to end on planted counts held in memory: S, G2/M and
//! resting cells binned apart, the model putting cells back in their bins,
//! a scored subset keeping other rows out, and the flags' defaults.

use super::*;
use nalgebra_sparse::CooMatrix;

/// Counts in memory, all rows returned per block.
struct Mem {
    x: CscMatrix<f32>,
    names: Vec<Box<str>>,
}

impl ColumnBlocks for Mem {
    fn n_cols(&self) -> usize {
        self.x.ncols()
    }
    fn row_names(&self) -> anyhow::Result<Vec<Box<str>>> {
        Ok(self.names.clone())
    }
    fn block(&self, lb: usize, ub: usize) -> anyhow::Result<CscMatrix<f32>> {
        let mut coo = CooMatrix::new(self.x.nrows(), ub - lb);
        for j in lb..ub {
            let col = self.x.col(j);
            for (&r, &v) in col.row_indices().iter().zip(col.values()) {
                coo.push(r, j - lb, v);
            }
        }
        Ok(CscMatrix::from(&coo))
    }
}

/// The Tirosh S and G2/M genes, 400 others and, with `peaks`, 200 more rows
/// of heavy counts (another modality's); 30 cells each of S, G2/M and
/// resting.
fn planted(peaks: bool) -> Mem {
    let cycle = crate::cell_cycle::CycleGenes::tirosh_2016();
    let s: Vec<&str> = cycle.symbols(crate::cell_cycle::CyclePhase::S).collect();
    let g: Vec<&str> = cycle.symbols(crate::cell_cycle::CyclePhase::G2M).collect();
    let mut names: Vec<Box<str>> = s.iter().chain(&g).map(|n| Box::from(*n)).collect();
    let n_cycle = names.len();
    names.extend((0..400).map(|i| format!("OTHER{i}").into_boxed_str()));
    let n_rna = names.len();
    if peaks {
        names.extend((0..200).map(|i| format!("chr1:{i}-{}", i + 1).into_boxed_str()));
    }
    let mut coo = CooMatrix::new(names.len(), 90);
    for c in 0..90 {
        for k in 0..400 {
            coo.push(n_cycle + k, c, 1.0 + ((k + c) % 8) as f32);
        }
        for r in 0..n_cycle {
            let own = match c / 30 {
                0 => r < s.len(),
                1 => r >= s.len(),
                _ => false,
            };
            coo.push(r, c, if own { 8.0 } else { 1.0 });
        }
        if peaks {
            // Peak counts that differ by cell, unrelated to the cycle.
            for p in 0..200 {
                coo.push(n_rna + p, c, 50.0 + ((p * 7 + c * 3) % 40) as f32);
            }
        }
    }
    Mem {
        x: CscMatrix::from(&coo),
        names,
    }
}

fn args(names: &[&str]) -> CovariateArgs {
    CovariateArgs {
        covariates: names.iter().map(|n| Box::from(*n)).collect(),
        ..CovariateArgs::default()
    }
}

/// A projection with no structure of its own.
fn flat_projection(n: usize) -> DMatrix<f32> {
    DMatrix::from_fn(4, n, |d, i| ((i * 7 + d * 13) % 11) as f32)
}

use nalgebra::DMatrix;

fn bin(l: &str) -> Option<&str> {
    l.split('|').nth(1)
}

/// The bin most of `cells`' binned cells have.
fn major(labels: &[Box<str>], cells: std::ops::Range<usize>) -> Option<&str> {
    let mut count: BTreeMap<&str, usize> = BTreeMap::new();
    for i in cells {
        if let Some(b) = bin(&labels[i]) {
            *count.entry(b).or_default() += 1;
        }
    }
    count.into_iter().max_by_key(|(_, c)| *c).map(|(b, _)| b)
}

use std::collections::BTreeMap;

#[test]
fn planted_states_take_bins_of_their_own_and_the_model_bins_them_again() -> anyhow::Result<()> {
    let data = planted(false);
    let batch: Vec<Box<str>> = vec!["b0".into(); 90];
    let m = model_covariates(
        &data,
        &batch,
        &flat_projection(90),
        &args(&["cycle"]),
        None,
        16,
    )?
    .expect("modelled");
    let (s, g, r) = (
        major(&m.labels, 0..30),
        major(&m.labels, 30..60),
        major(&m.labels, 60..90),
    );
    assert!(s.is_some() && g.is_some() && r.is_some());
    assert!(s != g && g != r && s != r, "{:?}", m.labels);
    assert_eq!(m.scores.ncols(), 2);
    // The recorded model, after a trip through a file, bins the cells the
    // fit binned the same way.
    let dir = tempfile::tempdir()?;
    let path = dir.path().join("m.json");
    m.model.write(&path)?;
    let back = CovariateModel::read(&path)?.expect("written");
    let again = back.labels(&data, &batch, None, 16)?;
    for (a, l) in again.iter().zip(&m.labels) {
        if bin(l).is_some() {
            assert_eq!(a, l);
        }
    }
    // Labels an earlier fit wrote are binned afresh.
    let refit = model_covariates(
        &data,
        &m.labels,
        &flat_projection(90),
        &args(&["cycle"]),
        None,
        16,
    )?
    .expect("modelled again");
    assert_eq!(refit.labels, m.labels);
    Ok(())
}

#[test]
fn a_scored_subset_keeps_other_rows_out_of_the_scores() -> anyhow::Result<()> {
    let rna = planted(false);
    let both = planted(true);
    let scored: Vec<usize> = (0..rna.names.len()).collect();
    let batch: Vec<Box<str>> = vec!["b0".into(); 90];
    let a = model_covariates(
        &rna,
        &batch,
        &flat_projection(90),
        &args(&["cycle"]),
        None,
        16,
    )?
    .expect("modelled");
    let b = model_covariates(
        &both,
        &batch,
        &flat_projection(90),
        &args(&["cycle"]),
        Some(&scored),
        16,
    )?
    .expect("modelled");
    // Scored on the RNA rows only, the peaks change nothing.
    assert!((a.scores.clone() - b.scores.clone()).abs().max() < 1e-4);
    assert_eq!(a.labels, b.labels);
    Ok(())
}

#[test]
fn nothing_is_modelled_without_genes_and_bad_flags_are_refused() -> anyhow::Result<()> {
    let data = Mem {
        x: planted(false).x,
        names: (0..497).map(|i| format!("G{i}").into_boxed_str()).collect(),
    };
    let batch: Vec<Box<str>> = vec!["b".into(); 90];
    let p = flat_projection(90);
    assert!(model_covariates(&data, &batch, &p, &args(&["all"]), None, 16)?.is_none());
    assert!(model_covariates(&data, &batch, &p, &args(&["none"]), None, 16)?.is_none());
    assert!(model_covariates(&data, &batch, &p, &args(&["nonsense"]), None, 16).is_err());
    assert!(model_covariates(&data, &batch, &p, &args(&["all", "cycle"]), None, 16).is_err());
    assert!(model_covariates(&data, &batch[..10], &p, &args(&["all"]), None, 16).is_err());
    Ok(())
}

#[test]
fn labels_with_parts_of_their_own_leave_all_be_and_refuse_a_list() -> anyhow::Result<()> {
    let data = planted(false);
    let batch: Vec<Box<str>> = vec!["donor1|lane2".into(); 90];
    let p = flat_projection(90);
    assert!(model_covariates(&data, &batch, &p, &args(&["all"]), None, 16)?.is_none());
    assert!(model_covariates(&data, &batch, &p, &args(&["cycle"]), None, 16).is_err());
    Ok(())
}

#[test]
fn new_cells_are_scored_on_the_rows_the_fit_scored() -> anyhow::Result<()> {
    let rna = planted(false);
    let both = planted(true);
    let batch: Vec<Box<str>> = vec!["b0".into(); 90];
    let m = model_covariates(
        &rna,
        &batch,
        &flat_projection(90),
        &args(&["cycle"]),
        None,
        16,
    )?
    .expect("modelled");
    assert_eq!(m.model.rows, rna.names);
    // Rows the fit never saw (another modality's) leave the bins as they were.
    let again = m.model.labels(&both, &batch, None, 16)?;
    for (a, l) in again.iter().zip(&m.labels) {
        if bin(l).is_some() {
            assert_eq!(a, l);
        }
    }
    Ok(())
}

#[test]
fn labels_compose_and_strip_one_covariate_axis() {
    assert_eq!(plain_batch("rep1").unwrap().as_ref(), "rep1");
    assert_eq!(plain_batch("rep1|cv3").unwrap().as_ref(), "rep1");
    assert!(plain_batch("rep1|donor").is_err());
    assert!(plain_batch("rep1|cv").is_err());
    assert!(plain_batch("a|b|cv0").is_err());
    assert_eq!(compose("rep1", Some(2)).unwrap().as_ref(), "rep1|cv2");
    assert_eq!(compose("rep1", None).unwrap().as_ref(), "rep1");
    assert_eq!(bin_of("cv12"), Some(12));
}

#[test]
fn an_old_run_replays_without_covariates_and_the_cli_defaults_to_all_and_icl() {
    let replayed: CovariateArgs = serde_json::from_str("{}").unwrap();
    assert!(!replayed.is_on().unwrap());
    #[derive(clap::Parser)]
    struct Cli {
        #[command(flatten)]
        c: CovariateArgs,
    }
    let cli = <Cli as clap::Parser>::parse_from(["x"]);
    assert!(cli.c.is_on().unwrap());
    assert_eq!(cli.c.covariates, [Box::from("all")]);
    assert_eq!(cli.c.covariate_criterion, Criterion::Icl);
    assert_eq!(cli.c.covariate_bins_max, 10);
}

#[test]
fn scoring_then_binning_is_the_one_call_and_bin_checks_its_cells() -> anyhow::Result<()> {
    let data = planted(false);
    let batch: Vec<Box<str>> = vec!["b0".into(); 90];
    let p = flat_projection(90);
    let one = model_covariates(&data, &batch, &p, &args(&["cycle"]), None, 16)?.expect("modelled");
    let scored = || score_covariates(&data, &args(&["cycle"]), None, 16);
    let s = scored()?.expect("scored");
    assert_eq!(s.scores.nrows(), 90);
    let two = s.bin(&batch, &p)?.expect("binned");
    assert_eq!(one.labels, two.labels);
    assert_eq!(one.model.rows, two.model.rows);
    assert!(scored()?.expect("scored").bin(&batch[..10], &p).is_err());
    assert!(scored()?
        .expect("scored")
        .bin(&batch, &flat_projection(10))
        .is_err());
    Ok(())
}

#[test]
fn a_named_covariate_without_its_genes_is_refused_and_all_scores_nothing() -> anyhow::Result<()> {
    let data = Mem {
        x: planted(false).x,
        names: (0..497).map(|i| format!("G{i}").into_boxed_str()).collect(),
    };
    assert!(score_covariates(&data, &args(&["all"]), None, 16)?.is_none());
    assert!(score_covariates(&data, &args(&["mito"]), None, 16).is_err());
    Ok(())
}

#[test]
fn stale_covariate_files_are_listed_and_cleared() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    let prefix = dir.path().join("run").to_string_lossy().into_owned();
    for p in [
        files::model(&prefix),
        files::cells(&prefix),
        files::batch(&prefix, 0),
        files::batch(&prefix, 1),
    ] {
        std::fs::write(p, "x")?;
    }
    assert_eq!(files::all(&prefix).len(), 4);
    assert_eq!(files::clear_stale(&prefix)?.len(), 4);
    assert!(files::all(&prefix).is_empty());
    assert!(files::is_batch(&files::batch(&prefix, 12)));
    assert!(files::is_batch("run.covariate.batch.0.txt"));
    assert!(!files::is_batch("run.batch.0.txt"));
    assert!(!files::is_batch(&files::model(&prefix)));
    Ok(())
}
