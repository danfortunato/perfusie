function [S, meta] = superedge_graph(G, superedges, opts)
%SUPEREDGE_GRAPH Build a coarse incidence graph from super-edge paths.
%   [S, META] = SUPEREDGE_GRAPH(G, SUPEREDGES) converts SUPEREDGES from
%   SMOOTH_GRAPH into a graph whose nodes are super-vertices and whose edges
%   are super-edges.
%
%   Super-vertices are:
%     - all degree ~= 2 nodes in G (anchors), and
%     - one representative node for each cycle super-edge.
%
%   Edge direction follows traversal in SUPEREDGES(k).nodes:
%     start = nodes(1), stop = nodes(end).
%
%   Name-value options:
%     opts.directed                 - return a digraph (default: true)
%     opts.include_isolated_anchors - include degree~=2 anchors even when
%                                     they are not part of any super-edge
%                                     (default: false)
%
%   Output edge variables in S:
%     superedge_id   - index into SUPEREDGES
%     orig_start     - original node index at traversal start
%     orig_end       - original node index at traversal end
%     is_cycle       - true for cycle super-edges
%     num_path_nodes - number of original nodes along the super-edge
%     path_nodes     - cell array of original node-index paths
%     path_length_xy - super-edge polyline length in x/y
%
%   Output node variables in S:
%     orig_id          - original node index in G
%     degree_original  - node degree in G
%     is_anchor        - degree ~= 2 in G

arguments
    G {mustBeA(G, 'graph')}
    superedges struct
    opts.directed (1, 1) logical = true
    opts.include_isolated_anchors (1, 1) logical = false
end

n = numnodes(G);
deg = degree(G);

if ~isempty(superedges) && ~isfield(superedges, 'nodes')
    error('superedge_graph:badSuperedges', ...
        'SUPEREDGES must contain field ''nodes''.');
end

has_coords = ismember('coords', G.Nodes.Properties.VariableNames);
if has_coords
    coords = G.Nodes.coords;
    if size(coords, 1) ~= n || size(coords, 2) < 2
        error('superedge_graph:badCoords', ...
            'G.Nodes.coords must be [numnodes(G) x >=2].');
    end
else
    coords = zeros(0, 2);
end

n_superedges = numel(superedges);
edge_orig_start = zeros(n_superedges, 1);
edge_orig_end = zeros(n_superedges, 1);
edge_is_cycle = false(n_superedges, 1);
edge_nodes = cell(n_superedges, 1);
edge_length = nan(n_superedges, 1);
cycle_reps = zeros(0, 1);

for k = 1:n_superedges
    ids = superedges(k).nodes;
    if ~(isnumeric(ids) && isvector(ids) && ~isempty(ids))
        error('superedge_graph:badNodes', ...
            'SUPEREDGES(%d).nodes must be a non-empty numeric vector.', k);
    end

    ids = ids(:);
    if any(ids < 1 | ids > n | fix(ids) ~= ids)
        error('superedge_graph:badNodeIndex', ...
            'SUPEREDGES(%d).nodes contains invalid node indices.', k);
    end

    is_cycle = false;
    if isfield(superedges, 'is_cycle')
        is_cycle = logical(superedges(k).is_cycle);
    end

    % Normalize cycle paths: keep exactly one copy of each node in traversal.
    if is_cycle && numel(ids) >= 2 && ids(1) == ids(end)
        ids = ids(1:end-1);
    end
    if isempty(ids)
        error('superedge_graph:emptyCyclePath', ...
            'SUPEREDGES(%d).nodes became empty after cycle normalization.', k);
    end

    edge_is_cycle(k) = is_cycle;
    edge_nodes{k} = ids;

    if is_cycle
        edge_orig_start(k) = ids(1);
        edge_orig_end(k) = ids(1);
        cycle_reps(end + 1, 1) = ids(1); %#ok<AGROW>
    else
        edge_orig_start(k) = ids(1);
        edge_orig_end(k) = ids(end);
    end

    if isfield(superedges, 'length')
        lenk = superedges(k).length;
        if isnumeric(lenk) && isscalar(lenk) && isfinite(lenk)
            edge_length(k) = lenk;
        end
    end
    if isnan(edge_length(k)) && has_coords
        edge_length(k) = polyline_length_xy(coords, ids, is_cycle);
    end
end

anchor_orig = find(deg ~= 2);
edge_vertex_orig = [edge_orig_start edge_orig_end].';
edge_vertex_orig = unique(edge_vertex_orig(:), 'stable');
% edge_vertex_orig = unique([edge_orig_start; edge_orig_end], 'stable');
if opts.include_isolated_anchors
    supervertex_orig = unique([edge_vertex_orig; cycle_reps; anchor_orig], 'stable');
else
    supervertex_orig = unique([edge_vertex_orig; cycle_reps], 'stable');
end

if isempty(supervertex_orig) && n > 0
    % Fallback for degenerate inputs (for example, manually supplied empty SUPEREDGES).
    supervertex_orig = 1;
end

n_superverts = numel(supervertex_orig);
orig_to_supervertex = zeros(n, 1);
orig_to_supervertex(supervertex_orig) = (1:n_superverts)';

edge_super_start = orig_to_supervertex(edge_orig_start);
edge_super_end   = orig_to_supervertex(edge_orig_end);

if any(edge_super_start == 0 | edge_super_end == 0)
    error('superedge_graph:missingSupervertex', ...
        'At least one super-edge endpoint did not map to a super-vertex.');
end

edge_order = (1:n_superedges)';

EdgeTable = table( ...
  edge_order, edge_orig_start, edge_orig_end, edge_is_cycle, ...
  cellfun(@numel, edge_nodes), edge_nodes, edge_length, ...
  'VariableNames', {'superedge_id','orig_start','orig_end', ...
                    'is_cycle','num_path_nodes','path_nodes','path_length_xy'});

if opts.directed
    S = digraph(edge_super_start, edge_super_end, EdgeTable, n_superverts);
else
    S = graph(edge_super_start, edge_super_end, EdgeTable, n_superverts);
end

if numedges(S) ~= n_superedges
    error('superedge_graph:edgeCountMismatch', ...
        ['Constructed graph has %d edges, but SUPEREDGES has %d entries. ' ...
         'This usually indicates unsupported graph simplification.'], ...
        numedges(S), n_superedges);
end

S.Nodes.orig_id = supervertex_orig;
S.Nodes.degree_original = deg(supervertex_orig);
S.Nodes.is_anchor = deg(supervertex_orig) ~= 2;

if has_coords
    S.Nodes.coords = coords(supervertex_orig, :);
end
if ismember('radius', G.Nodes.Properties.VariableNames)
    S.Nodes.radius = G.Nodes.radius(supervertex_orig, :);
end

% S.Edges.superedge_id = (1:n_superedges)';
% S.Edges.orig_start = edge_orig_start;
% S.Edges.orig_end = edge_orig_end;
% S.Edges.is_cycle = edge_is_cycle;
% S.Edges.num_path_nodes = cellfun(@numel, edge_nodes);
% S.Edges.path_nodes = edge_nodes;
% S.Edges.path_length_xy = edge_length;

meta = struct();
meta.orig_to_supervertex = orig_to_supervertex;
meta.supervertex_orig = supervertex_orig;
meta.superedge_endnodes_super = [edge_super_start edge_super_end];
meta.superedge_endnodes_orig = [edge_orig_start edge_orig_end];
meta.superedge_nodes_orig = edge_nodes;

end

function L = polyline_length_xy(coords, ids, is_cycle)
if numel(ids) < 2
    L = 0;
    return;
end

p = coords(ids, 1:2);
dp = diff(p, 1, 1);
L = sum(hypot(dp(:, 1), dp(:, 2)));

if is_cycle
    L = L + hypot(p(end, 1) - p(1, 1), p(end, 2) - p(1, 2));
end
end
