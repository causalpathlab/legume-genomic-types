//! Batch labels with more than one axis: a cell's batch and, say, its
//! cell-cycle bin, joined into one label (`rep1_wt|cc1`) that every tool
//! reading batch labels treats as one batch, and that a tool knowing the
//! axes splits back.
//!
//! One separator, [`AXIS_SEP`], for every tool, so a label written by one
//! (senna's cell-cycle bins) reads the same in another (lupin, a table of
//! per-batch effects). A part holding the separator is refused rather than
//! silently split wrong.

/// The separator between a label's axes.
pub const AXIS_SEP: char = '|';

/// `parts` joined into one label: `["rep1_wt", "cc1"]` → `rep1_wt|cc1`.
///
/// # Errors
/// No parts, an empty part, or a part holding [`AXIS_SEP`].
pub fn compose(parts: &[&str]) -> Result<String, String> {
    if parts.is_empty() {
        return Err("a batch label needs at least one part".into());
    }
    for p in parts {
        if p.is_empty() {
            return Err("a batch label's part is empty".into());
        }
        if p.contains(AXIS_SEP) {
            return Err(format!(
                "batch label part {p:?} holds {AXIS_SEP:?}, which joins a label's axes; rename it"
            ));
        }
    }
    Ok(parts.join(&AXIS_SEP.to_string()))
}

/// A label's axes, in order; a plain label is its one axis.
#[must_use]
pub fn split(label: &str) -> Vec<&str> {
    label.split(AXIS_SEP).collect()
}

/// The first axis of a label: the batch itself, before any axis joined to it.
#[must_use]
pub fn base(label: &str) -> &str {
    label.split(AXIS_SEP).next().unwrap_or(label)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labels_compose_and_split_back() {
        let l = compose(&["rep1_wt", "cc1"]).unwrap();
        assert_eq!(l, "rep1_wt|cc1");
        assert_eq!(split(&l), ["rep1_wt", "cc1"]);
        assert_eq!(base(&l), "rep1_wt");
        assert_eq!(split("rep1_wt"), ["rep1_wt"]);
        assert_eq!(base("rep1_wt"), "rep1_wt");
    }

    #[test]
    fn a_part_holding_the_separator_or_nothing_is_refused() {
        assert!(compose(&["a|b", "cc0"]).unwrap_err().contains("holds"));
        assert!(compose(&["", "cc0"]).is_err());
        assert!(compose(&[]).is_err());
    }
}
