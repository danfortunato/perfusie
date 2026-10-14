%% Smooth and plot the graph

close all
G = graphtools.gallery('L2_1_Tr9L');
smooth = 1.5;
fittol = 1e-4;

figure(1), clf, graphtools.plot_graph(G, root=true, vertex_number=true)

[G, superedges] = graphtools.smooth_graph(G, smooth);
figure(2), clf, graphtools.plot_graph(G, root=true, vertex_number=true)

cgrph = graphtools.graph_to_chunkgraph(G, superedges, tol=fittol);
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

drawnow
figure(1)
