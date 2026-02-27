function S = coarsen_graph(G, max_separation, relabel_opt, root)
%COARSEN_GRAPH MATLAB port of Python get_features.coarsen_graph.
%   S = COARSEN_GRAPH(G, MAX_SEPARATION, RELABEL_OPT, ROOT)
%   mirrors terminal-cells/data_analysis_code/get_features.py.
%
%   Inputs:
%     G              - undirected MATLAB graph with node variable `coords`
%     max_separation - segment length threshold (default: inf)
%     relabel_opt    - if true, reorder node IDs by distance from root (default: true)
%     root           - root node index in G (default: 1)
%
%   Output:
%     S              - coarsened graph. Node variable `orig_id` stores the
%                      corresponding node index in G.

if nargin < 2 || isempty(max_separation)
    max_separation = inf;
end
if nargin < 3 || isempty(relabel_opt)
    relabel_opt = true;
end
if nargin < 4 || isempty(root)
    root = 1;
end

if ~isa(G, 'graph')
    error('coarsen_graph:badInput', 'G must be a MATLAB graph.');
end
if ~ismember('coords', G.Nodes.Properties.VariableNames)
    error('coarsen_graph:missingCoords', 'G.Nodes.coords is required.');
end

n = numnodes(G);
if root < 1 || root > n || fix(root) ~= root
    error('coarsen_graph:badRoot', 'root must be an integer node index in [1, numnodes(G)].');
end
if ~(isnumeric(max_separation) && isscalar(max_separation) && ...
        (isinf(max_separation) || max_separation > 0))
    error('coarsen_graph:badMaxSeparation', ...
        'max_separation must be a positive scalar or inf.');
end

edges_g = edge_endnodes_numeric(G);
coords = G.Nodes.coords;
if size(coords, 1) ~= n || size(coords, 2) < 2
    error('coarsen_graph:badCoords', 'G.Nodes.coords must be [numnodes(G) x >=2].');
end

% Build adjacency list of G (node labels are original indices 1..n).
adj = build_adjacency(n, edges_g);

% First stage: fully simplified graph H (remove all degree-2 non-root nodes).
active = true(n, 1);
deg0 = degree_from_adj(adj, active);
special_count = nnz(deg0 ~= 2 | (1:n)' == root);

while nnz(active) > special_count
    deg = degree_from_adj(adj, active);
    idx = find(active & deg == 2 & (1:n)' ~= root, 1, 'first');
    if isempty(idx)
        break;
    end

    neibs = active_neighbors(adj{idx}, active);
    if numel(neibs) >= 2
        a = neibs(1);
        b = neibs(2);

        if ~ismember(b, adj{a})
            adj{a}(end + 1) = b;
        end
        if ~ismember(a, adj{b})
            adj{b}(end + 1) = a;
        end

        adj{a}(adj{a} == idx) = [];
        adj{b}(adj{b} == idx) = [];
    end

    old_neibs = adj{idx};
    for k = 1:numel(old_neibs)
        nb = old_neibs(k);
        adj{nb}(adj{nb} == idx) = [];
    end

    adj{idx} = [];
    active(idx) = false;
end

h_nodes = find(active);
h_edges = adjacency_to_edges(adj, active);

if isinf(max_separation)
    s_nodes = h_nodes;
    s_edges = h_edges;
else
    endpoint_merge_ratio = 0.5;
    new_nodes = zeros(0, 1);
    new_edges = zeros(0, 2);

    for i = 1:size(h_edges, 1)
        e0 = h_edges(i, 1);
        e1 = h_edges(i, 2);
        path = shortestpath(G, e0, e1);
        if isempty(path)
            continue;
        end

        shortened_path = path(1);
        n_anchor = e0;

        for j = 1:numel(path)
            m = path(j);
            c0 = coords(n_anchor, 1:2);
            c1 = coords(m, 1:2);
            d = hypot(c0(1) - c1(1), c0(2) - c1(2));

            if d > max_separation
                shortened_path(end + 1) = m; %#ok<AGROW>
                n_anchor = m;
            end
        end

        if shortened_path(end) ~= path(end)
            shortened_path(end + 1) = path(end);
        end

        shortened_path = prune_short_endpoint_segments( ...
            shortened_path, coords, max_separation, endpoint_merge_ratio);

        new_nodes = [new_nodes; shortened_path(:)]; %#ok<AGROW>
        if numel(shortened_path) > 1
            new_edges = [new_edges; ...
                [shortened_path(1:end-1)' shortened_path(2:end)']]; %#ok<AGROW>
        end
    end

    s_nodes = unique(new_nodes, 'stable');

    % Python behavior: S = G.subgraph(nodes).copy(); then add new_edges.
    in_sub = ismember(edges_g(:, 1), s_nodes) & ismember(edges_g(:, 2), s_nodes);
    sub_edges = edges_g(in_sub, :);
    s_edges = [sub_edges; new_edges];
    s_edges = unique_undirected_edges(s_edges);
end

S = build_graph_from_original_nodes(G, s_nodes, s_edges);

if relabel_opt && numnodes(S) > 0
    root_in_s = find(S.Nodes.orig_id == root, 1);
    if isempty(root_in_s)
        error('coarsen_graph:rootDropped', 'Root node %d was not present in coarsened graph.', root);
    end

    dist = distances(S, root_in_s)';
    path_len = dist + 1;
    tie = (1:numnodes(S))';
    [~, order] = sortrows([path_len tie], [1 2]);
    S = reordernodes(S, order);
end
end

function path = prune_short_endpoint_segments(path, coords, max_separation, ratio)
% Remove interior nodes that create very short end segments near endpoints.
if numel(path) < 3
    return;
end

min_sep = ratio * max_separation;

changed = true;
while changed && numel(path) >= 3
    changed = false;

    d_start = pair_dist_xy(coords, path(1), path(2));
    if d_start < min_sep && numel(path) >= 3
        path(2) = [];
        changed = true;
        if numel(path) < 3
            break;
        end
    end

    d_end = pair_dist_xy(coords, path(end - 1), path(end));
    if d_end < min_sep && numel(path) >= 3
        path(end - 1) = [];
        changed = true;
    end
end
end

function d = pair_dist_xy(coords, a, b)
c0 = coords(a, 1:2);
c1 = coords(b, 1:2);
d = hypot(c0(1) - c1(1), c0(2) - c1(2));
end

function edges = edge_endnodes_numeric(G)
end_nodes = G.Edges.EndNodes;
if isnumeric(end_nodes)
    edges = end_nodes;
elseif iscell(end_nodes) || isstring(end_nodes)
    if isstring(end_nodes)
        s = end_nodes(:, 1);
        t = end_nodes(:, 2);
    else
        s = string(end_nodes(:, 1));
        t = string(end_nodes(:, 2));
    end
    edges = [findnode(G, s), findnode(G, t)];

    if any(edges(:) < 1)
        error('coarsen_graph:badNames', 'Could not resolve named edge endpoints.');
    end
else
    error('coarsen_graph:badEdges', 'Unsupported EndNodes format.');
end
end

function adj = build_adjacency(n, edges)
adj = cell(n, 1);
for i = 1:size(edges, 1)
    u = edges(i, 1);
    v = edges(i, 2);
    adj{u}(end + 1) = v;
    adj{v}(end + 1) = u;
end
end

function d = degree_from_adj(adj, active)
n = numel(adj);
d = zeros(n, 1);
for i = 1:n
    if ~active(i)
        continue;
    end
    d(i) = nnz(active(adj{i}));
end
end

function nb = active_neighbors(neighbor_list, active)
if isempty(neighbor_list)
    nb = zeros(0, 1);
    return;
end
nb = neighbor_list(active(neighbor_list));
end

function edges = adjacency_to_edges(adj, active)
edges = zeros(0, 2);
for u = 1:numel(adj)
    if ~active(u)
        continue;
    end
    nb = adj{u};
    nb = nb(active(nb));
    nb = nb(nb > u);
    if ~isempty(nb)
        edges = [edges; [u * ones(numel(nb), 1) nb(:)]]; %#ok<AGROW>
    end
end
end

function e = unique_undirected_edges(e)
if isempty(e)
    return;
end
u = min(e(:, 1), e(:, 2));
v = max(e(:, 1), e(:, 2));
e = unique([u v], 'rows', 'stable');
end

function S = build_graph_from_original_nodes(G, orig_nodes, orig_edges)
orig_nodes = orig_nodes(:);
n = numel(orig_nodes);

if n == 0
    S = graph();
    S.Nodes.orig_id = zeros(0, 1);
    return;
end

map = zeros(numnodes(G), 1);
map(orig_nodes) = 1:n;

if isempty(orig_edges)
    s = zeros(0, 1);
    t = zeros(0, 1);
else
    s = map(orig_edges(:, 1));
    t = map(orig_edges(:, 2));
end

S = graph(s, t, [], n);
S.Nodes.orig_id = orig_nodes;

node_vars = G.Nodes.Properties.VariableNames;
for i = 1:numel(node_vars)
    var_name = node_vars{i};
    if strcmp(var_name, 'orig_id')
        continue;
    end
    S.Nodes.(var_name) = G.Nodes.(var_name)(orig_nodes, :);
end
end
