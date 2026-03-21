function cgrph = fit(file, opts)

arguments
    file
    opts.smooth = 1
    opts.order = 16
    opts.tol = 1e-6
end

G = graphtools.import_graph(file);
[G, superedges] = graphtools.smooth_graph(G, opts.smooth);
cgrph = graphtools.graph_to_chunkgraph(G, superedges, tol=opts.tol, order=opts.order);

end
