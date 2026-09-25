# legume-genomic-types

Shared genomic types and parsers for legume tools (`senna`, `chickpea`,
`faba`, [`mung-cnv`](https://github.com/causalpathlab/mung-cnv), …):

- GFF/GTF parsing
- BED intervals
- SAM/BAM barcode and strand helpers
- Transcript / exon models

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
