function G = gallery(varargin)

if ( nargin == 1 )
    n = varargin{1};
    files = dir('terminal-cells/data/traces_L*/*.traces');
    file = fullfile(files(n).folder, files(n).name);
elseif ( nargin == 2 )
    [folder, name] = deal(varargin{:});
    file = fullfile('terminal-cells', 'data', folder, [name '.traces']);
end

G = graphtools.import_graph(file);

end
