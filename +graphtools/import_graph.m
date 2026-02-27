function G = import_graph(filename)
%IMPORT_GRAPH   Read a terminal-cell .traces file into a MATLAB graph.
%   G = IMPORT_GRAPH(FILENAME) mirrors the Python routine
%   trace_file_to_G in terminal-cells/data_analysis_code/convert_trace_to_network.py.
%
%   Output graph fields:
%     G.Nodes.coords  - [N x 3] node coordinates (x,y,z), with y-axis flipped
%     G.Nodes.radius  - [N x 1] node radius
%     G.Edges.length  - [E x 1] edge length in projected 2D (x,y)

arguments
    filename {mustBeTextScalar}
end
filename = char(filename);

if ~isfile(filename)
    error('trace_file_to_G:fileNotFound', 'File not found: %s', filename);
end

xml_file = filename;
tmp_dir = '';
cleanup_obj = onCleanup(@() []);

if is_gzip_file(filename)
    tmp_dir = tempname;
    mkdir(tmp_dir);
    cleanup_obj = onCleanup(@() cleanup_tmp_dir(tmp_dir));
    out_files = gunzip(filename, tmp_dir);
    if isempty(out_files)
        error('trace_file_to_G:gunzipFailed', 'Could not decompress: %s', filename);
    end
    xml_file = out_files{1};
end

doc = xmlread(xml_file);
root = doc.getDocumentElement();
children = root.getChildNodes();

coords = zeros(256, 3);
radius = zeros(256, 1);
edge_s = zeros(512, 1);
edge_t = zeros(512, 1);
node_count = 1;
edge_count = 0;

for i = 0:(children.getLength() - 1)
    path_node = children.item(i);

    if path_node.getNodeType() ~= org.w3c.dom.Node.ELEMENT_NODE
        continue;
    end
    if ~strcmp(char(path_node.getNodeName()), 'path')
        continue;
    end

    [path_keys, path_vals] = ordered_attributes(path_node);
    fit_idx = find(strcmp(path_keys, 'usefitted'), 1);
    if isempty(fit_idx) || ~strcmpi(path_vals{fit_idx}, 'false')
        continue;
    end

    point_nodes = path_node.getElementsByTagName('point');
    first_pt = true;

    for j = 0:(point_nodes.getLength() - 1)
        point_node = point_nodes.item(j);
        [~, point_vals] = ordered_attributes(point_node);

        if numel(point_vals) < 6
            error('trace_file_to_G:pointFormat', ...
                'Unexpected point attribute format in %s.', filename);
        end


        x = str2double(point_node.getAttribute('xd'));
        y = str2double(point_node.getAttribute('yd'));
        z = str2double(point_node.getAttribute('zd'));
        r = str2double(point_node.getAttribute('x'));

        ensure_node_capacity();
        coords(node_count, :) = [x y z];
        radius(node_count) = r;

        if node_count > 1
            if first_pt
                neighbor = closest_node(coords(1:(node_count - 1), :), [x, y, z]);
                first_pt = false;
            else
                neighbor = node_count - 1;
            end

            ensure_edge_capacity();
            edge_count = edge_count + 1;
            edge_s(edge_count) = node_count;
            edge_t(edge_count) = neighbor;
        end

        node_count = node_count + 1;
    end
end

n_nodes = node_count - 1;
coords = coords(1:n_nodes, :);
radius = radius(1:n_nodes, :);

if n_nodes > 0
    max_y = max(coords(:, 2));
    coords(:, 2) = max_y - coords(:, 2);
end

if edge_count > 0
    edge_s = edge_s(1:edge_count);
    edge_t = edge_t(1:edge_count);
    G = graph(edge_s, edge_t);
else
    G = graph();
end

if numnodes(G) < n_nodes
    G = addnode(G, n_nodes - numnodes(G));
end

G.Nodes.coords = coords;
G.Nodes.radius = radius;

if edge_count > 0
    c1 = coords(edge_s, 1:2);
    c2 = coords(edge_t, 1:2);
    edge_len = hypot(c1(:, 1) - c2(:, 1), c1(:, 2) - c2(:, 2));
else
    edge_len = zeros(0, 1);
end
G.Edges.length = edge_len;

clear cleanup_obj; %#ok<CLSCR>

    function ensure_node_capacity()
        if node_count > size(coords, 1)
            new_n = max(2 * size(coords, 1), node_count);
            coords(end + 1:new_n, :) = 0;
            radius(end + 1:new_n, :) = 0;
        end
    end

    function ensure_edge_capacity()
        next_edge = edge_count + 1;
        if next_edge > size(edge_s, 1)
            new_n = max(2 * size(edge_s, 1), next_edge);
            edge_s(end + 1:new_n, :) = 0;
            edge_t(end + 1:new_n, :) = 0;
        end
    end
end

function idx = closest_node(existing_coords, query_coord)
% Return index of closest existing node (first occurrence on ties).
sq_dist = sum((existing_coords - query_coord) .^ 2, 2);
[~, idx] = min(sq_dist);
end

function [keys, vals] = ordered_attributes(node)
attrs = node.getAttributes();
n = attrs.getLength();
keys = cell(1, n);
vals = cell(1, n);
for k = 1:n
    attr = attrs.item(k - 1);
    keys{k} = char(attr.getName());
    vals{k} = char(attr.getValue());
end
end

function tf = is_gzip_file(filename)
fid = fopen(filename, 'r');
if fid < 0
    error('trace_file_to_G:openFailed', 'Could not open file: %s', filename);
end
cleaner = onCleanup(@() fclose(fid));
sig = fread(fid, 2, '*uint8');
tf = numel(sig) == 2 && sig(1) == 31 && sig(2) == 139;
end

function cleanup_tmp_dir(tmp_dir)
if isfolder(tmp_dir)
    rmdir(tmp_dir, 's');
end
end
