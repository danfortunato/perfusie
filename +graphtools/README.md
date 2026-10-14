# graphtools

Tools for turning traced terminal-cell networks into chunkgraphs for
`perfusionop`. Several are MATLAB ports of Python routines in
[terminal-cells](https://github.com/tgavrilchenko/terminal-cells)
(`data_analysis_code/`).

- `gallery` loads a cell from the package's `data/` folder by name,
  e.g. `graphtools.gallery('L2_1_Tr9L')`.
- `import_graph` reads a `.traces` file into a MATLAB `graph` with node
  coordinates and radii (port of `trace_file_to_G`).
- `smooth_graph` splits a graph into super-edges (chains between branch
  points and tips) and smooths each one with a Gaussian along arc length.
- `superedge_graph` builds the coarse graph whose edges are the
  super-edges.
- `graph_to_chunkgraph` fits each super-edge with a chunker and assembles a
  chunkgraph.
- `fit` does all of the above in one call: smooth, then fit a chunkgraph.
- `coarsen_graph` removes degree-2 nodes, optionally keeping a maximum node
  spacing (port of `get_features.coarsen_graph`).
- `rasterize_graph` draws a graph's edges onto a logical image.
- `plot_graph` plots a graph at its spatial coordinates.
- `digraph_to_graph` converts a digraph to an undirected graph.

Typical use:

```matlab
G = graphtools.gallery('L2_1_Tr9L');
cgrph = graphtools.fit(G, smooth=1, tol=1e-3);
```
