function G = gallery(name)
%GALLERY   Load a terminal-cell network from the bundled data.
%   G = GALLERY(NAME) imports the cell NAME from the package's data/ folder
%   as a MATLAB graph, e.g. GALLERY('L2_1_Tr9L'). Cell names have the form
%   L<stage>_<larva>_Tr9<side>; see data/README.md. Each cell's .pdf image
%   in data/ shows its geometry.
%
%   G = GALLERY(N) imports the N-th cell in alphabetical order.
%
%   See also GRAPHTOOLS.IMPORT_GRAPH, GRAPHTOOLS.FIT.

datadir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data');

if ( isnumeric(name) )
    files = dir(fullfile(datadir, '*.traces'));
    file = fullfile(datadir, files(name).name);
else
    file = fullfile(datadir, [char(name) '.traces']);
    if ( ~isfile(file) )
        error('graphtools:gallery:notFound', ...
            'No cell named ''%s'' in %s.', name, datadir);
    end
end

G = graphtools.import_graph(file);

end
