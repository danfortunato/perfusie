function cgrph = graph_to_chunkgraph(G, superedges, opts)
%GRAPH_TO_CHUNKGRAPH   Fit a chunkgraph to the super-edges of a graph.
%   CGRPH = GRAPH_TO_CHUNKGRAPH(G, SUPEREDGES) fits each super-edge of G
%   (from GRAPHTOOLS.SMOOTH_GRAPH) with an open chunker through its node
%   coordinates G.Nodes.coords(:,1:2), and assembles them into a chunkgraph.
%   The chunkgraph's vertices are the super-vertices of G, numbered as in
%   GRAPHTOOLS.SUPEREDGE_GRAPH, and edge k is super-edge k, oriented from
%   SUPEREDGES(k).nodes(1) to SUPEREDGES(k).nodes(end).
%
%   CGRPH = GRAPH_TO_CHUNKGRAPH(G, SUPEREDGES, Name=Value) sets options:
%     tol   - tolerance for the chunker fit (default 1e-6).
%     order - Legendre order of each chunk (default 16).
%     nover - number of times to oversample (uniformly refine) each fitted
%             chunker (default 0).
%
%   See also GRAPHTOOLS.FIT, CHUNKERFIT, CHUNKGRAPH.

arguments
    G
    superedges
    opts.tol = 1e-6
    opts.order = 16
    opts.nover = 0
end

chnkrs = cell(length(superedges), 1);
chopts = [];
chopts.ifclosed = false;
chopts.cparams = [];
chopts.cparams.eps = opts.tol;
chopts.cparams.nover = opts.nover;
chopts.pref = [];
chopts.pref.k = opts.order;
for k = 1:length(superedges)
    nodes = G.Nodes.coords(superedges(k).nodes, 1:2);
    chnkrs{k} = chunkerfit(nodes.', chopts);
end

Gd = graphtools.superedge_graph(G, superedges, directed=true);
verts = Gd.Nodes.coords(:,1:2).';
edges = sortrows(Gd.Edges, 'superedge_id');
edgesendverts = edges.EndNodes.';
cgrph = chunkgraph(verts, edgesendverts, chnkrs);

end
