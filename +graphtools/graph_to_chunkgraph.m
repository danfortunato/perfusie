function cgrph = graph_to_chunkgraph(G, superedges, opts)

arguments
    G
    superedges
    opts.tol = 1e-6
    opts.order = 16
end

chnkrs = cell(length(superedges), 1);
chopts = [];
chopts.ifclosed = false;
chopts.cparams = [];
chopts.cparams.eps = opts.tol;
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
