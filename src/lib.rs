//! Genomic data structures and parsers
//!
//! This library provides reusable genomic data structures including:
//! - GFF/GTF parsing
//! - SAM/BAM utilities
//! - BED file format
//! - Genomic positions
//! - Variant ids, genome builds and a streaming VCF reader
//! - The E2G enhancer–gene release (parquet)

pub mod bed;
pub mod coordinates;
pub mod e2g;
pub mod gff;
pub mod plink;
pub mod positions;
pub mod sam;
pub mod transcript;
pub mod variant;
pub mod vcf;
