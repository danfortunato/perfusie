% Solve the perfusion model on a traced terminal cell, check O2
% conservation, and plot the tissue partial pressure.

%%% Graph parameters
name = 'L2_1_Tr9L';   % Terminal cell (see data/)
smooth = 1;
fittol = 1e-3;
n = 16;
nover = 1;
outer = true;

%%% Physical parameters
root   = 1;     % Indices of trachea roots
alpha  = 250;   % Scaled channel area (1/xi)
lambda = 0;     % Inverse permeability
phi    = 0.1;   % Yukawa parameter
eta    = 0.03;  % Henry's constant

%%% Solver parameters
method   = 'direct';
schur    = true;
fmm      = false;
precond  = abs(lambda) < 0.01;
showcond = false;
tol      = 1e-6;

% Import the graph, smooth it, and discretize it
t = tic;
cgrph = graphtools.fit(graphtools.gallery(name), smooth=smooth, tol=fittol, order=n, nover=nover);
fprintf('Time to process graph:      %gs\n', toc(t));

t = tic;
op = perfusionop(cgrph, root=root, outer=outer, alpha=alpha, lambda=lambda, phi=phi, eta=eta);
[u, vars] = op.solve(method=method, schur=schur, fmm=fmm, precond=precond, tol=tol, showcond=showcond);
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
