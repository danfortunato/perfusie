function [BW, xg, yg] = rasterize_graph(G, dx)
%RASTERIZE_GRAPH Rasterize a geometric graph onto a Cartesian grid.
%   BW = RASTERIZE_GRAPH(G, DX) returns a logical matrix BW whose true
%   entries approximate graph edges in G using XY coordinates from
%   G.Nodes.coords(:,1:2) on a grid with spacing DX.
%
%   [BW, XG, YG] also returns grid coordinate vectors so plotting with
%   imagesc(XG, YG, BW); axis xy; gives Cartesian orientation.

arguments
    G {mustBeA(G, 'graph')}
    dx (1, 1) double {mustBePositive}
end

if ~ismember('coords', G.Nodes.Properties.VariableNames)
    error('rasterize_graph:missingCoords', 'G.Nodes.coords is required.');
end

xy = G.Nodes.coords;
if size(xy, 2) < 2
    error('rasterize_graph:badCoords', 'G.Nodes.coords must have at least 2 columns.');
end
xy = xy(:, 1:2);
n = numnodes(G);
if size(xy, 1) ~= n
    error('rasterize_graph:badCoords', 'G.Nodes.coords row count must match numnodes(G).');
end

if n == 0
    BW = false(0, 0);
    xg = zeros(0, 1);
    yg = zeros(0, 1);
    return;
end

xmin = min(xy(:, 1)) - 0.5 * dx;
xmax = max(xy(:, 1)) + 0.5 * dx;
ymin = min(xy(:, 2)) - 0.5 * dx;
ymax = max(xy(:, 2)) + 0.5 * dx;

nx = max(1, ceil((xmax - xmin) / dx) + 1);
ny = max(1, ceil((ymax - ymin) / dx) + 1);
xg = xmin + (0:(nx - 1)) * dx;
yg = ymin + (0:(ny - 1)) * dx;

BW = false(ny, nx);
end_nodes = numeric_endnodes(G);

% Rasterize each edge by sampling at <= dx/2 along the segment.
for k = 1:size(end_nodes, 1)
    i = end_nodes(k, 1);
    j = end_nodes(k, 2);
    p0 = xy(i, :);
    p1 = xy(j, :);
    seg_len = hypot(p1(1) - p0(1), p1(2) - p0(2));
    n_samp = max(2, ceil(seg_len / (0.5 * dx)) + 1);
    t = linspace(0, 1, n_samp);
    xs = p0(1) + t * (p1(1) - p0(1));
    ys = p0(2) + t * (p1(2) - p0(2));
    [rows, cols] = xy_to_sub(xs, ys, xmin, ymin, dx, ny, nx);
    BW(sub2ind([ny, nx], rows, cols)) = true;
end

% Ensure all vertices are represented in the raster.
[rows_v, cols_v] = xy_to_sub(xy(:, 1), xy(:, 2), xmin, ymin, dx, ny, nx);
BW(sub2ind([ny, nx], rows_v, cols_v)) = true;
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
    error('rasterize_graph:badEndNodes', 'Unsupported edge endpoint type.');
end
edges = [findnode(G, s), findnode(G, t)];
if any(edges(:) < 1)
    error('rasterize_graph:badEndNodes', 'Could not resolve edge endpoints.');
end
end

function [rows, cols] = xy_to_sub(x, y, xmin, ymin, dx, ny, nx)
cols = round((x - xmin) ./ dx) + 1;
rows = round((y - ymin) ./ dx) + 1;
cols = min(max(cols, 1), nx);
rows = min(max(rows, 1), ny);
rows = rows(:);
cols = cols(:);
end
