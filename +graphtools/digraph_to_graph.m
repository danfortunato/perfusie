function G = digraph_to_graph(D)
%DIGRAPH_TO_GRAPH   Convert a digraph to an undirected graph.
%   G = DIGRAPH_TO_GRAPH(D) returns the undirected graph with an edge
%   wherever D has an edge in either direction, keeping the node table of
%   D. Edge properties of D are not kept.

A = adjacency(D);
G = graph(A+A.', D.Nodes);

end
