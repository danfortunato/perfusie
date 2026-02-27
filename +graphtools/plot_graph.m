function plot_graph(G, opts)

arguments
    G
    opts.root          = false
    opts.vertex_number = false
    opts.edge_number   = false
end

xy = G.Nodes.coords(:,1:2);
plot(G, xdata=xy(:,1), ydata=xy(:,2), nodecolor='k', edgecolor='k', nodelabel={})

hold on

if ( opts.root )
    plot(xy(1,1), xy(1,2), 'ro', 'LineWidth', 1.5)
end

if ( opts.vertex_number )
    is_directed = isa(G, 'digraph');
    for k = 1:numnodes(G)
        if is_directed
            d = numel(unique([predecessors(G, k); successors(G, k)]));
        else
            d = numel(neighbors(G, k));
        end
        if d ~= 2
            text(xy(k,1), xy(k,2), [' ' num2str(k)], color='b', fontsize=10)
        end
    end
end

if ( opts.edge_number )
    for k = 1:numedges(G)
        v = G.Edges.EndNodes(k,:);
        mid = 0.5*(xy(v(1),:) + xy(v(2),:));
        if ismember('superedge_id', fieldnames(G.Edges))
            id = G.Edges.superedge_id(k);
        else
            id = k;
        end
        text(mid(1), mid(2), [' ' num2str(id)], color='r', fontsize=10)
    end
end

hold off

axis equal tight
xlabel('x')
ylabel('y')

end
