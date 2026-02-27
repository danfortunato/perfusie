%%% Graph parameters
file = fullfile('terminal-cells', 'data', 'traces_L2', '1_Tr9L.traces');
smooth = 1;

%%% Physical parameters
xi = 0.001;   % Inverse channel width scaled by diffusion ratio
alpha = 0.1;  % Inverse permeability
phi = 0.1;    % Yukawa parameter
root = 1;     % Indices of trachea roots

%%% Solver parameters
method = 'gmres';
schur  = true;
fmm    = true;

% Import the graph, smooth it, and discretize it
t = tic;
G = graphtools.import_graph(file);
[G, superedges] = graphtools.smooth_graph(G, smooth);
cgrph = graphtools.graph_to_chunkgraph(G, superedges);
fprintf('Time to process graph:     %gs\n', toc(t));

t = tic;
op = perfusionop(cgrph, xi=xi, alpha=alpha, phi=phi, root=root);
u = op.solve(method=method, schur=schur, fmm=fmm);
fprintf('Time to solve:             %gs\n', toc(t));

t = tic;
figure(1), clf
plot_perfusion(cgrph, u)
fprintf('Time to evaluate and plot: %gs\n', toc(t));
colormap turbo
shg
