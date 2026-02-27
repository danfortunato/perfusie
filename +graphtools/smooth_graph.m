function [G_smooth, superedges] = smooth_graph(G, sigma, opts)
%SMOOTH_GRAPH   Apply 1D Gaussian smoothing along super-edges.
%   [G_SMOOTH, SUPEREDGES] = SMOOTH_GRAPH(G, SIGMA, OPTS)
%   smooths geometry using G.Nodes.coords(:,1:2).
%
%   A super-edge is a maximal chain where interior nodes have degree 2 and
%   endpoints have degree ~= 2 (or a closed degree-2 cycle).
%
%   Required:
%     G     - MATLAB graph with node variable `coords`
%     sigma - Gaussian standard deviation in the same units as coordinates

%   Optional (name-value via opts):
%     opts.truncate        - kernel radius multiplier (default: 3)
%     opts.fix_endpoints   - keep super-edge endpoints fixed (default: true)
%     opts.taper_endpoints - taper smoothing near endpoints (default: false)
%     opts.redistribute    - re-space nodes along smoothed curve (default: true)
%     opts.equal_spacing   - final pass: equal arc-length node spacing on
%                            each super-edge (default: true)

arguments
    G {mustBeA(G, 'graph')}
    sigma (1, 1) double {mustBePositive}
    opts.truncate (1, 1) double {mustBePositive} = 3
    opts.fix_endpoints (1, 1) logical = true
    opts.taper_endpoints (1, 1) logical = false
    opts.redistribute (1, 1) logical = true
    opts.equal_spacing (1, 1) logical = true
end

if ~ismember('coords', G.Nodes.Properties.VariableNames)
    error('smooth_superedges_gaussian:missingCoords', 'G.Nodes.coords is required.');
end

coords = G.Nodes.coords;
n = numnodes(G);
if size(coords, 1) ~= n || size(coords, 2) < 2
    error('smooth_superedges_gaussian:badCoords', ...
        'G.Nodes.coords must be [numnodes(G) x >=2].');
end

xy = coords(:, 1:2);
edges = numeric_endnodes(G);
deg = degree(G);

superedges = extract_superedges(edges, deg, n, xy);
xy_new = xy;

for k = 1:numel(superedges)
    ids = superedges(k).nodes;
    if numel(ids) < 2
        continue;
    end

    if superedges(k).is_cycle
        xy_seg_new = smooth_cycle_segment(xy(ids, :), sigma, opts.truncate, ...
            opts.redistribute, opts.equal_spacing);
        write_mask = true(size(ids));
    else
        xy_seg_new = smooth_open_segment(xy(ids, :), sigma, ...
            opts.truncate, opts.fix_endpoints, opts.taper_endpoints, ...
            opts.redistribute, opts.equal_spacing);

        write_mask = true(size(ids));
        if opts.fix_endpoints && numel(ids) >= 2
            write_mask([1, end]) = false;
        end
    end

    if any(write_mask)
        xy_new(ids(write_mask), :) = xy_seg_new(write_mask, :);
    end
end

coords_new = coords;
coords_new(:, 1:2) = xy_new;

G_smooth = G;
G_smooth.Nodes.coords = coords_new;

if ismember('length', G_smooth.Edges.Properties.VariableNames) && ~isempty(edges)
    c1 = xy_new(edges(:, 1), :);
    c2 = xy_new(edges(:, 2), :);
    G_smooth.Edges.length = hypot(c1(:, 1) - c2(:, 1), c1(:, 2) - c2(:, 2));
end
end

function segments = extract_superedges(edges, deg, n, xy)
[adj_nodes, adj_edges] = adjacency_with_edge_ids(edges, n);
visited = false(size(edges, 1), 1);
segments = struct('nodes', {}, 'is_cycle', {}, 'length', {});

anchors = find(deg ~= 2);

for ia = 1:numel(anchors)
    a = anchors(ia);
    for k = 1:numel(adj_nodes{a})
        e0 = adj_edges{a}(k);
        if visited(e0)
            continue;
        end

        nxt = adj_nodes{a}(k);
        nodes = [a; nxt];
        visited(e0) = true;

        prev = a;
        curr = nxt;
        while deg(curr) == 2
            idx = find(adj_nodes{curr} ~= prev, 1, 'first');
            if isempty(idx)
                break;
            end
            nn = adj_nodes{curr}(idx);
            ee = adj_edges{curr}(idx);
            if visited(ee)
                break;
            end
            nodes(end + 1, 1) = nn; %#ok<AGROW>
            visited(ee) = true;
            prev = curr;
            curr = nn;
        end

        seg.nodes = nodes;
        seg.is_cycle = false;
        seg.length = polyline_length(nodes, false);
        segments(end + 1) = seg; %#ok<AGROW>
    end
end

% Remaining unvisited edges belong to degree-2 cycles.
while true
    e0 = find(~visited, 1, 'first');
    if isempty(e0)
        break;
    end

    u = edges(e0, 1);
    v = edges(e0, 2);
    nodes = [u; v];
    visited(e0) = true;

    prev = u;
    curr = v;
    while true
        idx = find(adj_nodes{curr} ~= prev, 1, 'first');
        if isempty(idx)
            break;
        end
        nn = adj_nodes{curr}(idx);
        ee = adj_edges{curr}(idx);

        if visited(ee)
            break;
        end

        nodes(end + 1, 1) = nn; %#ok<AGROW>
        visited(ee) = true;
        prev = curr;
        curr = nn;

        if curr == u
            break;
        end
    end

    if numel(nodes) >= 2 && nodes(end) == nodes(1)
        nodes(end) = [];
    end

    seg.nodes = nodes;
    seg.is_cycle = true;
    seg.length = polyline_length(nodes, true);
    segments(end + 1) = seg; %#ok<AGROW>
end

    function L = polyline_length(node_ids, close_loop)
        if numel(node_ids) < 2
            L = 0;
            return;
        end
        p = xy(node_ids, :);
        dp = diff(p, 1, 1);
        L = sum(hypot(dp(:, 1), dp(:, 2)));
        if close_loop
            L = L + hypot(p(end, 1) - p(1, 1), p(end, 2) - p(1, 2));
        end
    end
end

function [adj_nodes, adj_edges] = adjacency_with_edge_ids(edges, n)
adj_nodes = cell(n, 1);
adj_edges = cell(n, 1);
for e = 1:size(edges, 1)
    u = edges(e, 1);
    v = edges(e, 2);
    adj_nodes{u}(end + 1) = v;
    adj_edges{u}(end + 1) = e;
    adj_nodes{v}(end + 1) = u;
    adj_edges{v}(end + 1) = e;
end
end

function P_new = smooth_open_segment(P, sigma, truncate, fix_endpoints, taper_endpoints, redistribute, equal_spacing)
if size(P, 1) < 3
    P_new = P;
    return;
end

s = [0; cumsum(hypot(diff(P(:, 1)), diff(P(:, 2))))];
L = s(end);
if L <= 0
    P_new = P;
    return;
end

% Tangent smoothing preserves segment lengths and avoids center collapse.
P_new = smooth_open_by_tangent(P, sigma, truncate);

if redistribute %#ok<NASGU>
    % No-op: tangent reconstruction already preserves node spacing.
end

if taper_endpoints
    t = s / L;
    taper = sin(pi * t) .^ 2;
    P_new = P + taper .* (P_new - P);
end

if fix_endpoints
    P_new = enforce_open_endpoints(P_new, P);
end

if equal_spacing
    P_new = redistribute_open_equal(P_new, size(P, 1));
    if fix_endpoints
        P_new([1, end], :) = P([1, end], :);
    end
end
end

function P_new = smooth_cycle_segment(P, sigma, truncate, redistribute, equal_spacing)
if size(P, 1) < 3
    P_new = P;
    return;
end

P_new = smooth_cycle_by_tangent(P, sigma, truncate);

if redistribute %#ok<NASGU>
    % No-op: tangent reconstruction already preserves node spacing.
end

if equal_spacing
    P_new = redistribute_cycle_equal(P_new, size(P, 1));
end
end

function P_sm = smooth_open_by_tangent(P, sigma, truncate)
V = diff(P, 1, 1);
lens = hypot(V(:, 1), V(:, 2));
if ~any(lens > 0)
    P_sm = P;
    return;
end

s_nodes = [0; cumsum(lens)];
L = s_nodes(end);
seg_pos = s_nodes(1:end-1) + 0.5 * lens;

T = unit_rows(V);
T_sm = smooth_tangents_open(T, seg_pos, sigma, truncate);
P_sm = reconstruct_from_tangents_open(P(1, :), lens, T_sm);

% Distribute endpoint correction along arclength so ends match exactly.
delta = P(end, :) - P_sm(end, :);
tau = s_nodes / max(L, eps);
P_sm = P_sm + tau .* delta;
end

function P_sm = smooth_cycle_by_tangent(P, sigma, truncate)
V = [diff(P, 1, 1); P(1, :) - P(end, :)];
lens = hypot(V(:, 1), V(:, 2));
if ~any(lens > 0)
    P_sm = P;
    return;
end

s_nodes = [0; cumsum(lens)];
L = s_nodes(end);
seg_pos = s_nodes(1:end-1) + 0.5 * lens;

T = unit_rows(V);
T_sm = smooth_tangents_cycle(T, seg_pos, L, sigma, truncate);
Pc = reconstruct_from_tangents_cycle(P(1, :), lens, T_sm);

% Enforce closure with distributed correction.
err = Pc(end, :) - Pc(1, :);
tau = s_nodes / max(L, eps);
Pc = Pc - tau .* err;
Pc(end, :) = Pc(1, :);

P_sm = Pc(1:end-1, :);
end

function T_sm = smooth_tangents_open(T, seg_pos, sigma, truncate)
n = size(T, 1);
T_sm = zeros(n, 2);
cutoff = truncate * sigma;

for i = 1:n
    d = abs(seg_pos - seg_pos(i));
    use = d <= cutoff;
    if ~any(use)
        T_sm(i, :) = T(i, :);
        continue;
    end
    w = exp(-0.5 * (d(use) / sigma) .^ 2);
    u = (w' * T(use, :)) / sum(w);
    nu = hypot(u(1), u(2));
    if nu <= 1e-14
        T_sm(i, :) = T(i, :);
    else
        T_sm(i, :) = u / nu;
    end
end
end

function T_sm = smooth_tangents_cycle(T, seg_pos, L, sigma, truncate)
n = size(T, 1);
T_sm = zeros(n, 2);
cutoff = truncate * sigma;

for i = 1:n
    d = abs(seg_pos - seg_pos(i));
    d = min(d, L - d);
    use = d <= cutoff;
    if ~any(use)
        T_sm(i, :) = T(i, :);
        continue;
    end
    w = exp(-0.5 * (d(use) / sigma) .^ 2);
    u = (w' * T(use, :)) / sum(w);
    nu = hypot(u(1), u(2));
    if nu <= 1e-14
        T_sm(i, :) = T(i, :);
    else
        T_sm(i, :) = u / nu;
    end
end
end

function P_rec = reconstruct_from_tangents_open(p0, lens, T)
nseg = numel(lens);
P_rec = zeros(nseg + 1, 2);
P_rec(1, :) = p0;
for i = 1:nseg
    P_rec(i + 1, :) = P_rec(i, :) + lens(i) * T(i, :);
end
end

function Pc = reconstruct_from_tangents_cycle(p0, lens, T)
nseg = numel(lens);
Pc = zeros(nseg + 1, 2);
Pc(1, :) = p0;
for i = 1:nseg
    Pc(i + 1, :) = Pc(i, :) + lens(i) * T(i, :);
end
end

function P_out = enforce_open_endpoints(P_in, P_ref)
n = size(P_in, 1);
if n < 2
    P_out = P_in;
    return;
end

s = [0; cumsum(hypot(diff(P_in(:, 1)), diff(P_in(:, 2))))];
L = s(end);
if L <= 0
    P_out = P_in;
    P_out([1, end], :) = P_ref([1, end], :);
    return;
end

delta0 = P_ref(1, :) - P_in(1, :);
delta1 = P_ref(end, :) - P_in(end, :);
tau = s / L;
P_out = P_in + (1 - tau) .* delta0 + tau .* delta1;
P_out([1, end], :) = P_ref([1, end], :);
end

function P_eq = redistribute_open_equal(P, npt)
if npt <= 1
    P_eq = P;
    return;
end

s = [0; cumsum(hypot(diff(P(:, 1)), diff(P(:, 2))))];
L = s(end);
if L <= 0
    P_eq = P;
    return;
end

[su, Pu] = unique_curve_samples(s, P);
if numel(su) < 2
    P_eq = P;
    return;
end

targ = linspace(0, L, npt)';
x = interp1(su, Pu(:, 1), targ, 'linear', 'extrap');
y = interp1(su, Pu(:, 2), targ, 'linear', 'extrap');
P_eq = [x, y];
end

function P_eq = redistribute_cycle_equal(P, npt)
if npt <= 2
    P_eq = P;
    return;
end

Pc = [P; P(1, :)];
s = [0; cumsum(hypot(diff(Pc(:, 1)), diff(Pc(:, 2))))];
L = s(end);
if L <= 0
    P_eq = P;
    return;
end

[su, Pu] = unique_curve_samples(s, Pc);
if numel(su) < 2
    P_eq = P;
    return;
end
if su(end) < L
    su(end + 1, 1) = L;
    Pu(end + 1, :) = Pu(1, :);
elseif norm(Pu(end, :) - Pu(1, :)) > 0
    su(end + 1, 1) = L;
    Pu(end + 1, :) = Pu(1, :);
end

targ = (0:npt-1)' * (L / npt);
x = interp1(su, Pu(:, 1), targ, 'linear', 'extrap');
y = interp1(su, Pu(:, 2), targ, 'linear', 'extrap');
P_eq = [x, y];
end

function [su, Pu] = unique_curve_samples(s, P)
tol = 1e-14;
keep = [true; diff(s) > tol];
su = s(keep);
Pu = P(keep, :);

if su(end) < s(end)
    su(end + 1, 1) = s(end);
    Pu(end + 1, :) = P(end, :);
end
end

function U = unit_rows(V)
U = V;
nr = hypot(V(:, 1), V(:, 2));
good = nr > 1e-14;
U(good, :) = V(good, :) ./ nr(good);

if any(~good)
    idx_good = find(good, 1, 'first');
    if isempty(idx_good)
        U(~good, :) = repmat([1, 0], nnz(~good), 1);
    else
        fallback = U(idx_good, :);
        U(~good, :) = repmat(fallback, nnz(~good), 1);
    end
end
end

function edges = numeric_endnodes(G)
end_nodes = G.Edges.EndNodes;
if isnumeric(end_nodes)
    edges = end_nodes;
    return;
end

if iscell(end_nodes)
    s = string(end_nodes(:, 1));
    t = string(end_nodes(:, 2));
elseif isstring(end_nodes)
    s = end_nodes(:, 1);
    t = end_nodes(:, 2);
else
    error('smooth_superedges_gaussian:badEndNodes', ...
        'Unsupported edge endpoint type.');
end

edges = [findnode(G, s), findnode(G, t)];
if any(edges(:) < 1)
    error('smooth_superedges_gaussian:badEndNodes', ...
        'Could not resolve named edge endpoints.');
end
end
