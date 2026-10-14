% Self-convergence test with adaptive (dyadic) refinement into the vertices
% of the graph. Vertices and branchpoints are where the singularities of the
% solution live, so resolving them with progressively finer panels (rather
% than uniformly refining the whole graph as in test_convergence.m) should
% give much faster convergence.
%
% At each refinement level we split the panel touching each vertex (i.e.
% the first and last panel of every edge) and let the chunker level
% restriction balance adjacent panels. The arclength of the vertex panels
% therefore halves at every level, while the smooth interior is left alone.
% With a high-order discretization we expect roughly spectral convergence
% in the vertex panel size h = h0 * 2^-k.

%%% Graph parameters
file = fullfile('terminal-cells', 'data', 'traces_L1', '1_Tr9L.traces');
smooth = 1;
fittol = 1e-3;
order  = 16;
outer  = false;

%%% Physical parameters
root  = 1;     % Indices of trachea roots
xi    = 0.004; % Inverse channel width scaled by diffusion ratio
alpha = 0.1;   % Inverse permeability
phi   = 0.1;   % Yukawa parameter
P     = 0.03;  % Henry's constant

%%% Solver parameters
method  = 'gmres';
schur   = true;
fmm     = true;
precond = abs(alpha) < 0.01;
tol     = 1e-10;

% Import the graph, smooth it, and discretize it
G = graphtools.import_graph(file);
[G, superedges] = graphtools.smooth_graph(G, smooth);
cgrph0 = graphtools.graph_to_chunkgraph(G, superedges, tol=fittol, order=order);

% Target points: one far from the network (smooth region), one close to a
% vertex (where the vertex singularity dominates the error)
xtarg_far   = 90;
ytarg_far   = 15;
xtarg_close = 98.5;
ytarg_close = 16.1;

figure(1), clf
plot(cgrph0)
hold on
plot(xtarg_close, ytarg_close, 'r.', markersize=15)
plot(xtarg_far,   ytarg_far,   'b.', markersize=15)
hold off
title('Initial graph and target points')
shg

%% Sweep over refinement levels. At each level we split the two panels of
%  every edge that touch a vertex (level restriction balances neighbours).

nref = 8;
utarg_far   = zeros(nref+1, 1);
utarg_close = zeros(nref+1, 1);
npts        = zeros(nref+1, 1);
hvals       = zeros(nref+1, 1);

cgrph = cgrph0;
for k = 1:nref+1
    npts(k)  = cgrph.npt;
    hvals(k) = min_vertex_arc(cgrph);
    fprintf('Level %d: h = %.3g, %d unknowns\n', k-1, hvals(k), npts(k));

    op = perfusionop(cgrph, root=root, outer=outer, xi=xi, alpha=alpha, phi=phi, P=P);
    [u, ~] = op.solve(method=method, schur=schur, fmm=fmm, precond=precond, tol=tol);

    utarg_far(k)   = u.bulk(xtarg_far,   ytarg_far);
    utarg_close(k) = u.bulk(xtarg_close, ytarg_close);

    if ( k < nref+1 )
        cgrph = refine_into_vertices(cgrph);
    end
end

%% Plot the chunkgraph at the finest refinement level

figure(2), clf
plot(cgrph)
hold on
plot(xtarg_close, ytarg_close, 'r.', markersize=15)
plot(xtarg_far,   ytarg_far,   'b.', markersize=15)
hold off
title(sprintf('Finest refinement: h = %.3g, npt = %d', hvals(end), npts(end)))
shg

%% Convergence plot. The reference solution is the finest panelization.

err_far   = abs(utarg_far(1:end-1)   - utarg_far(end))   / abs(utarg_far(end));
err_close = abs(utarg_close(1:end-1) - utarg_close(end)) / abs(utarg_close(end));
h         = hvals(1:end-1);

figure(3), clf
loglog(1./h, err_close, '-ro', linewidth=1, markerfacecolor='r')
hold on
loglog(1./h, err_far,   '-bs', linewidth=1, markerfacecolor='b')
hold off
grid on
axis tight
legend( ...
    sprintf('Error at (%g,%g) (close to graph)', xtarg_close, ytarg_close), ...
    sprintf('Error at (%g,%g) (far from graph)', xtarg_far,   ytarg_far), ...
    interpreter='latex', fontsize=12, location='southwest')
xlabel('$1/h$ (inverse vertex-panel arclength)', interpreter='latex', fontsize=18)
ylabel('$\|u - u_\mathrm{ref}\| / \|u_\mathrm{ref}\|$', interpreter='latex', fontsize=18)
title(sprintf('Adaptive vertex refinement (order = %d)', order), fontsize=14)
shg

%% Convergence vs. number of unknowns

figure(4), clf
loglog(npts(1:end-1), err_close, '-ro', linewidth=1, markerfacecolor='r')
hold on
loglog(npts(1:end-1), err_far,   '-bs', linewidth=1, markerfacecolor='b')
hold off
grid on
axis tight
legend( ...
    sprintf('Error at (%g,%g) (close to graph)', xtarg_close, ytarg_close), ...
    sprintf('Error at (%g,%g) (far from graph)', xtarg_far,   ytarg_far), ...
    interpreter='latex', fontsize=12, location='southwest')
xlabel('number of unknowns', interpreter='latex', fontsize=18)
ylabel('$\|u - u_\mathrm{ref}\| / \|u_\mathrm{ref}\|$', interpreter='latex', fontsize=18)
title(sprintf('Adaptive vertex refinement (order = %d)', order), fontsize=14)
shg

%% Save plots to PNG

exportgraphics(figure(1), 'convergence_adaptive_graph_initial.png', resolution=200)
exportgraphics(figure(2), 'convergence_adaptive_graph_finest.png',  resolution=200)
exportgraphics(figure(3), 'convergence_adaptive_vs_h.png',          resolution=200)
exportgraphics(figure(4), 'convergence_adaptive_vs_npt.png',        resolution=200)

%% Local helpers

% Split the first and last panel of every edge. The chunkgraph refine
% routine will then enforce level restriction, so adjacent panels balance
% smoothly. Each call halves the arclength of every vertex panel.
function cg = refine_into_vertices(cg)
    splitchunks = cell(1, length(cg.echnks));
    for ie = 1:length(cg.echnks)
        nch_i = cg.echnks(ie).nch;
        if nch_i == 1
            splitchunks{ie} = 1;
        else
            splitchunks{ie} = [1, nch_i];
        end
    end
    cg = refine(cg, struct('splitchunks', {splitchunks}));
end

% Minimum arclength of any panel touching a vertex.
function h = min_vertex_arc(cg)
    arcs = [];
    for ie = 1:length(cg.echnks)
        w = sum(cg.echnks(ie).wts, 1);
        arcs = [arcs, w(1), w(end)]; %#ok<AGROW>
    end
    h = min(arcs);
end
