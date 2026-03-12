%%% Graph parameters
file = fullfile('terminal-cells', 'data', 'traces_L2', '1_Tr9L.traces');
smooth = 1;
outer = true;

%%% Physical parameters
root  = 1;     % Indices of trachea roots
xi    = 0.004; % Inverse channel width scaled by diffusion ratio
alpha = 0;     % Inverse permeability
phi   = 0.1;   % Yukawa parameter
P     = 0.03;  % Henry's constant

%%% Solver parameters
method  = 'gmres';
schur   = true;
fmm     = true;
precond = abs(alpha) < 0.01;
tol     = 1e-6;

% Import the graph, smooth it, and discretize it
t = tic;
G = graphtools.import_graph(file);
[G, superedges] = graphtools.smooth_graph(G, smooth);
cgrph = graphtools.graph_to_chunkgraph(G, superedges);
fprintf('Time to process graph:      %gs\n', toc(t));

t = tic;
op = perfusionop(cgrph, root=root, outer=outer, xi=xi, alpha=alpha, phi=phi, P=P);
[u, vars] = op.solve(method=method, schur=schur, fmm=fmm, precond=precond, tol=tol);
fprintf('Time to solve:              %gs\n', toc(t));

t = tic;
op.check_conservation(u, vars, ngrid=400);
fprintf('Time to check conservation: %gs\n', toc(t));

figure(1), clf
t = tic;
plot_perfusion(op, u, partialpressure=true, forcesmooth=true, nplot=800);
fprintf('Time to evaluate and plot:  %gs\n', toc(t));
colormap turbo
shg
