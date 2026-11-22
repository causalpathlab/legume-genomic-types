use super::*;

fn names(v: &[&str]) -> Vec<Box<str>> {
    v.iter().map(|s| Box::from(*s)).collect()
}

#[test]
fn gene_index_resolves_hgnc_renames_both_ways() {
    // A matrix on the CURRENT HGNC release; a marker panel on an older one. Without the
    // alias tier every one of these silently drops out of the panel.
    let dict = names(&[
        "H4C3", "H2BC21", "H1-0", "MARCHF2", "SEPTIN7", "MTARC1", "CD8A",
    ]);
    let idx = GeneIndex::build(&dict);
    assert_eq!(idx.match_gene("HIST1H4C"), Some(0));
    assert_eq!(idx.match_gene("HIST2H2BE"), Some(1));
    assert_eq!(idx.match_gene("H1F0"), Some(2));
    assert_eq!(idx.match_gene("MARCH2"), Some(3), "rule-based MARCH family");
    assert_eq!(idx.match_gene("SEPT7"), Some(4), "rule-based SEPT family");
    assert_eq!(
        idx.match_gene("MARC1"),
        Some(5),
        "MARC1 is MTARC1, not a MARCH-family member"
    );

    // ...and the reverse: an old-release matrix, a current panel.
    let old = names(&["HIST1H4C", "MARCH2", "SEPT7", "CD8A"]);
    let idx_old = GeneIndex::build(&old);
    assert_eq!(idx_old.match_gene("H4C3"), Some(0));
    assert_eq!(idx_old.match_gene("MARCHF2"), Some(1));
    assert_eq!(idx_old.match_gene("SEPTIN7"), Some(2));

    // The alias tier resolves through the faba `{gene}/count/{track}` row keys and the
    // `ENSG…_SYMBOL` form too, since it runs on the extracted symbol.
    let faba = names(&["ENSG00000197061_H4C3/count/spliced", "CD8A/count/spliced"]);
    let idx_faba = GeneIndex::build(&faba);
    assert_eq!(idx_faba.match_gene("HIST1H4C"), Some(0));

    // A gene with no alias must not be coerced into one.
    assert_eq!(idx.match_gene("ZZZ9"), None);
}

#[test]
fn gene_index_tiers_and_idf() {
    let dict = names(&["ENSG00000153563_CD8A", "MS4A1", "A_FOO", "FOO"]);
    let idx = GeneIndex::build(&dict);

    // exact (case-insensitive)
    assert_eq!(idx.match_gene("ms4a1"), Some(1));
    // symbol = last `_`-segment
    assert_eq!(idx.match_gene("CD8A"), Some(0));
    // exact preferred over an earlier-indexed suffix match (the refinement)
    assert_eq!(idx.match_gene("FOO"), Some(3));
    // unmatched
    assert_eq!(idx.match_gene("ZZZ9"), None);

    // bare ENSG query resolves to the `ENSG…_SYMBOL` row (and the combined
    // form resolves too) — HGNC / ENSG / ENSG_HGNC reconcile both ways.
    assert_eq!(idx.match_gene("ENSG00000153563"), Some(0));
    assert_eq!(idx.match_gene("ENSG00000153563_CD8A"), Some(0));

    // dict keyed by bare ENSG: a combined query still resolves via the
    // leading-id tier.
    let dict2 = names(&["ENSG00000153563", "MS4A1"]);
    let idx2 = GeneIndex::build(&dict2);
    assert_eq!(idx2.match_gene("ENSG00000153563_CD8A"), Some(0));
    assert_eq!(idx2.match_gene("ensg00000153563"), Some(0));

    // IDF: ubiquitous (df == C) → 0; exclusive → ln(C)
    assert_eq!(idf_weight(4, 4), 0.0);
    assert!(idf_weight(4, 1) > 1.38 && idf_weight(4, 1) < 1.39);
}

#[test]
fn loci_match_loci_strictly_by_key() {
    let dict = names(&["X:0-100", "chr1:0-100", "GENE1", "ENSG000_GENE2"]);
    let idx = GeneIndex::build(&dict);
    assert_eq!(idx.match_gene("chrX:0-100"), Some(0));
    assert_eq!(idx.match_gene("CHR1:0-100"), Some(1));
    assert_eq!(idx.match_gene("x:0-100"), None, "chromosome case is kept");
    assert_eq!(idx.match_gene("X:0-101"), None);
    assert_eq!(idx.match_gene("GENE2"), Some(3));
    // A gene query never lands on a locus row through the fallback tiers.
    assert_eq!(idx.match_gene("X"), None);
    assert_eq!(idx.match_gene("chr1"), None);
    // A locus query still reaches a row whose `/`-core is that locus.
    let idx = GeneIndex::build(&names(&["GENE1", "chr2:5-9/count/spliced"]));
    assert_eq!(idx.match_gene("chr2:5-9"), Some(1));
    assert_eq!(idx.match_gene("chr2:5-9/count/spliced"), Some(1));
}

#[test]
fn locus_queries_prefer_the_exact_row_and_reach_composite_rows() {
    let idx = GeneIndex::build(&names(&["chr2:5-9/count/spliced", "chr2:5-9"]));
    assert_eq!(idx.match_gene("chr2:5-9"), Some(1), "exact name first");
    assert_eq!(
        idx.match_gene("2:5-9"),
        Some(1),
        "a whole-locus row wins its key"
    );
    assert_eq!(idx.match_gene("chr2:5-9/count/spliced"), Some(0));
    assert_eq!(
        idx.match_gene("2:5-9/count/spliced"),
        Some(0),
        "chr prefix free"
    );
    let idx = GeneIndex::build(&names(&["chrX:0-100", "X:0-100"]));
    assert_eq!(idx.match_gene("X:0-100"), Some(1));
    assert_eq!(idx.match_gene("chrX:0-100"), Some(0));
    // An id_name composite peak row is reached by its locus.
    let idx = GeneIndex::build(&names(&["GENE9", "chr1:100-200_GENE1"]));
    assert_eq!(idx.match_gene("chr1:100-200"), Some(1));
    assert_eq!(idx.match_gene("x:0-100/count/spliced"), None);
}

/// The fallback tier finds the row a scan of the whole vocabulary finds.
#[test]
fn the_fallback_tier_agrees_with_a_full_scan() {
    let parts = ["a", "b", "ab", "c", "", "ensg0001", "x1", "gene"];
    // Names of one to four parts, thinned so that many queries have no row,
    // then a few with an aux suffix or a locus.
    let mut dict: Vec<String> = Vec::new();
    for n in 1..=4usize {
        for code in 0..parts.len().pow(n as u32) {
            let mut c = code;
            let segs: Vec<&str> = (0..n)
                .map(|_| {
                    let p = parts[c % parts.len()];
                    c /= parts.len();
                    p
                })
                .collect();
            if code % 7 < 2 {
                dict.push(segs.join("_"));
            }
        }
    }
    dict.push("A_B/count/spliced".into());
    dict.push("chr1:1-2".into());
    dict.push("chr1:1-2_A".into());
    let dict: Vec<Box<str>> = dict.into_iter().map(String::into_boxed_str).collect();
    let idx = GeneIndex::build(&dict);

    let scan = |gl: &str| {
        let suffix = format!("_{gl}");
        let prefix = format!("{gl}_");
        let middle = format!("_{gl}_");
        idx.lowered
            .iter()
            .position(|t| t.ends_with(&suffix) || t.starts_with(&prefix) || t.contains(&middle))
    };
    let mut queries: Vec<String> = vec![
        String::new(),
        "_".into(),
        "__".into(),
        "zzz".into(),
        "a_".into(),
        "_a".into(),
        "b/count/spliced".into(),
        "chr1:1-2".into(),
        "spliced".into(),
    ];
    for a in parts {
        queries.push(a.to_string());
        for b in parts {
            queries.push(format!("{a}_{b}"));
            for c in parts {
                queries.push(format!("{a}_{b}_{c}"));
            }
        }
    }
    let mut n_hits = 0;
    for q in &queries {
        let want = scan(q);
        assert_eq!(idx.match_segments(q), want, "query `{q}`");
        n_hits += usize::from(want.is_some());
    }
    assert!(n_hits > 50 && n_hits < queries.len(), "{n_hits} hits");
}
