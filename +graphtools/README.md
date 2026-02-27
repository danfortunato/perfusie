# Terminal-Cell Trace Reader (MATLAB)

- `graph_from_traces.m`
  - MATLAB port of Python `trace_file_to_G` from
    `terminal-cells/data_analysis_code/convert_trace_to_network.py`.
  - Reads gzipped XML `.traces` files and returns a `graph` with:
    - `G.Nodes.coords` (`[x y z]`, with y-axis flipped),
    - `G.Nodes.radius`,
    - `G.Edges.length` (2D projected edge length in `x,y`).

- `demo_read_trace_file.m`
  - Loads one sample file and plots the tree.

- `coarsen_graph.m`
  - MATLAB port of Python `get_features.coarsen_graph`.
  - Removes degree-2 nodes (except root), optionally enforces a
    `max_separation`, and optionally relabels by root distance.

- `rasterize_graph.m`
  - Rasterizes graph edges from `G.Nodes.coords(:,1:2)` to a logical grid
    with spacing `dx`.

- `smooth_superedges_gaussian.m`
  - Finds super-edges (degree-2 chains) and applies Gaussian smoothing
    along arc length per super-edge.

- `superedge_graph.m`
  - Builds a coarse superedge-incidence graph/digraph from `superedges`.
  - Each edge corresponds to one super-edge path, directed by path traversal.

## Usage

```matlab
run('matlab/setup.m');
G = graph_from_traces('terminal-cells/data/traces_L1/10_Tr9L.traces');
```
