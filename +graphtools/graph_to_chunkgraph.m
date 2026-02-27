function cgrph = graph_to_chunkgraph(G, superedges)

chnkrs = cell(length(superedges), 1);
opts = [];
opts.ifclosed = false;
for k = 1:length(superedges)
    nodes = G.Nodes.coords(superedges(k).nodes, 1:2);
    chnkrs{k} = chunkerfit(nodes.', opts);
end

Gd = graphtools.superedge_graph(G, superedges, directed=true);
verts = Gd.Nodes.coords(:,1:2).';
edges = sortrows(Gd.Edges, 'superedge_id');
edgesendverts = edges.EndNodes.';
cgrph = chunkgraph(verts, edgesendverts, chnkrs);

end
