function G = digraph_to_graph(D)

A = adjacency(D);
G = graph(A+A.', D.Nodes);

end
