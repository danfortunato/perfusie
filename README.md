# perfusie

A perfusion solver for 1D networks in 2D tissue using
[chunkIE](https://github.com/fastalgorithms/chunkie).

`perfusionop` solves a coupled PDE–ODE model of oxygen delivery through a
network of air-filled channels (such as an insect tracheal terminal cell)
into a thin slab of tissue: a modified Helmholtz equation in the tissue,
coupled through a permeable membrane to a 1D diffusion equation on each
network edge. The tissue concentration is represented as a single-layer
potential on the network, so only the network (and an optional outer
boundary) is discretized.

## Installation

perfusie is a [mip](https://mip.sh) package (mip 1.2.0 or later). Its
dependencies, chunkie (with fmm2d and flam) and minboundsuite, are installed
automatically from the `mip-org/core` channel:

```matlab
mip install -e /path/to/perfusie
mip load perfusie
```

Use an editable install (`-e`) so that edits to the source take effect
immediately.

## Quick start

```matlab
G = graphtools.gallery('L2_1_Tr9L');          % A traced terminal cell
cgrph = graphtools.fit(G, smooth=1, tol=1e-3);  % Smooth chunkgraph fit
op = perfusionop(cgrph, xi=0.004, phi=0.1, outer=true);
[u, vars] = op.solve();
plot_perfusion(op, u)
```

See `help perfusionop` for the model, parameters, and solver options.

## Contents

- `perfusionop.m` — the solver.
- `plot_perfusion.m` — plot a solution in the tissue and network.
- `+graphtools/` — load, smooth, and fit traced networks
  (see `+graphtools/README.md`).
- `data/` — 162 traced terminal cells, each with a PDF image, from
  Gavrilchenko et al. (see `data/README.md`).
- `examples/` — demos and convergence studies. Load them with
  `mip load perfusie --with examples`.
- `tests/` — unit tests.

## Tests

```matlab
mip test perfusie
```

When MATLAB is run with `-batch`, it may crash on exit after tests that use
the FMM. This is a known issue with gfortran-built MEX files and does not
affect the results.
