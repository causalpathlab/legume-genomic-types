# legume-genomic-types

Shared genomic types and parsers for legume tools (`senna`, `chickpea`,
`faba`, [`mung-cnv`](https://github.com/causalpathlab/mung-cnv), …):

- GFF/GTF parsing
- BED intervals
- SAM/BAM barcode and strand helpers
- Transcript / exon models
- Gene-name matching across conventions (`name_matching::GeneIndex`)
- Batch labels with more than one axis (`batch_label`): `rep1_wt|cc1`,
  one separator for every tool, composed and split back
- Cell-cycle marker genes (`cell_cycle`): the S and G2/M sets of
  Tirosh et al. 2016, *Science* 352:189 (doi:10.1126/science.aad0501),
  bundled as `data/cell_cycle_tirosh2016.tsv` with Seurat's 2019 HGNC
  renames; a table of your own in the same form replaces them

The Rust library crate is named `genomic_data`:

```toml
genomic-data = { version = "0.4.1", package = "legume-genomic-types" }
```

```sh
cargo add legume-genomic-types
# then: use genomic_data::...
```

## License

MIT
