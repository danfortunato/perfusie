# Terminal-cell networks

Traced *Drosophila* tracheal terminal cells (dorsal Tr9 cells) from

> T. Gavrilchenko, A. G. Simpkins, T. Simpson, L. A. Barrett, P. Hansen,
> S. Y. Shvartsman, and J. Schottenfeld-Roames, "The Drosophila tracheal
> terminal cell as a model for branching morphogenesis."

copied from <https://github.com/tgavrilchenko/terminal-cells> (commit
`98152d5`). That repository also holds the micrographs, masks, and analysis
code, which are not included here.

## Files

Each cell is named `L<stage>_<larva>_Tr9<side>`:

- `<stage>` is the larval stage, 1–3. L1 cells are the simplest, L3 the most
  branched.
- `<larva>` numbers the larva within that stage.
- `<side>` is `L` or `R` for the left or right cell.

For each cell there are two files:

- `.traces` is the traced network (gzipped SNT XML), read by
  `graphtools.import_graph`. Load one with `graphtools.gallery('L2_1_Tr9L')`.
- `.pdf` is an image of the filled cell, for choosing a geometry.

## Changes from the source repository

- Traces come from `data/traces_L<stage>/<larva>_Tr9<side>.traces`, renamed
  with an `L<stage>_` prefix so all stages fit in one folder.
- PDFs come from `data/all_individual_fills/`, with their original names,
  except that `L1_29_Tr9{L,R}.pdf` are renamed to `L1_25_Tr9{L,R}.pdf` to
  match the traces they show. `L1_30_Tr9{L,R}.pdf` match no trace and are
  omitted.
- `L1_16_Tr9{L,R}` have no individual PDF. Both cells appear in
  `data/fills_L1/16_Tr9.pdf` in the source repository.
