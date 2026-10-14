function cgrph = fit(G, opts)
%FIT   Fit a smooth chunkgraph to a traced network.
%   CGRPH = FIT(G) smooths the graph G along its super-edges and fits each
%   super-edge with a high-order chunker, returning a chunkgraph whose
%   edges are the super-edges of G. G may be a MATLAB graph (for example,
%   from GRAPHTOOLS.GALLERY) or the name of a .traces file.
%
%   CGRPH = FIT(G, Name=Value) sets options:
%     smooth - Gaussian smoothing width, in the units of the coordinates
%              (default 1). See GRAPHTOOLS.SMOOTH_GRAPH.
%     order  - Legendre order of each chunk (default 16).
%     tol    - tolerance for the chunker fit (default 1e-6).
%     nover  - number of times to oversample (uniformly refine) the fitted
%              chunkers (default 0).
%
%   See also GRAPHTOOLS.SMOOTH_GRAPH, GRAPHTOOLS.GRAPH_TO_CHUNKGRAPH.

arguments
    G
    opts.smooth = 1
    opts.order = 16
    opts.tol = 1e-6
    opts.nover = 0
end

if ( ~isa(G, 'graph') )
    G = graphtools.import_graph(G);
end
[G, superedges] = graphtools.smooth_graph(G, opts.smooth);
cgrph = graphtools.graph_to_chunkgraph(G, superedges, tol=opts.tol, order=opts.order, nover=opts.nover);

end
