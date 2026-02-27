%% Smooth and plot the graph

close all
G = graphtools.gallery('traces_L1', '10_Tr9L');
figure(1), clf, graphtools.plot_graph(G, root=true, vertex_number=true)

sigma = 1;
[G, superedges] = graphtools.smooth_graph(G, sigma);
figure(2), clf, graphtools.plot_graph(G, root=true, vertex_number=true)

cgrph = graphtools.graph_to_chunkgraph(G, superedges);

ax = axis;
figure(3)
plot(cgrph, '-o')
xlabel('x')
ylabel('y')
box on
axis equal
axis(ax)
shg

figure(4), clf
Gd = graphtools.superedge_graph(G, superedges, directed=true);
graphtools.plot_graph(Gd, root=true, vertex_number=true, edge_number=true)

%% Check if the endpoints of the chunkers match the original graph vertices

err = 0;
for k_vert = 1:length(cgrph.echnks)
    chnkr = cgrph.echnks(k_vert);
    endpts = chunkends(chnkr, [1 chnkr.nch]);
    vs1 = endpts(:,[1 4]);
    vs2 = G.Nodes.coords(superedges(k_vert).nodes([1 end]), 1:2).';
    err = err + norm(vs1 - vs2);
end
fprintf('Error = %g\n', err);

