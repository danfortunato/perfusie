% Test the convergence rate for a simple perfusion problem by looking at the
% value of u.bulk(x,y) at a fixed (x,y) under uniform refinement of the
% graph's panelization. We observe O(h^2) convergence due to the presence of
% endpoint/branchpoint singularities. (In practice, adaptive refinement should
% be performed at these locations, not uniform refinement.)

%%% Graph parameters
file = fullfile('terminal-cells', 'data', 'traces_L1', '1_Tr9L.traces');
smooth = 1;
outer = false;

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
cgrph = graphtools.graph_to_chunkgraph(G, superedges, tol=1e-3);

xtarg_far = 90;
ytarg_far = 15;
xtarg_close = 98.5;
ytarg_close = 16.1;

figure(1), clf
plot(cgrph)
hold on
plot(xtarg_close, ytarg_close, 'r.', markersize=10)
hold off
shg

%% Perform the sweep

nref = 9;
utarg_far = zeros(nref+1, 1);
utarg_close = zeros(nref+1, 1);
for k = 1:nref+1
    fprintf('Refinement level %d: %g unknowns\n', k-1, cgrph.npt);
    op = perfusionop(cgrph, root=root, outer=outer, xi=xi, alpha=alpha, phi=phi, P=P);
    [u, vars] = op.solve(method=method, schur=schur, fmm=fmm, precond=precond, tol=tol);
    utarg_far(k)   = u.bulk(xtarg_far,   ytarg_far);
    utarg_close(k) = u.bulk(xtarg_close, ytarg_close);
    if ( k < nref+1 )
        % Uniformly refine the graph once
        cgrph = refine(cgrph, struct(nover=1));
    end
end

%% Plot convergence of the error, taking the reference solution to be the
%  solution on the finest panelization

err_far   = abs(utarg_far(1:end-1)   - utarg_far(end))   / abs(utarg_far(end));
err_close = abs(utarg_close(1:end-1) - utarg_close(end)) / abs(utarg_close(end));
h = 2.^-(0:nref-1);

figure(2), clf
loglog(1./h, err_close, '-bo', linewidth=1)
hold on
loglog(1./h, 1e-6*h.^2, 'k--', linewidth=1)
hold off
grid on
axis tight
decades_equal(gca, 2)
legend( ...
    sprintf('Error at (%g,%g)', xtarg_close, ytarg_close), ...
    '$\mathcal{O}(h^2)$', ...
interpreter='latex', fontsize=16)
xlabel('$1/h$', interpreter='latex', fontsize=18)
ylabel('$\frac{\|u - u_\mathrm{ref}\|}{\|u_{ref}\|}$', interpreter='latex', fontsize=18)
shg

function decades_equal(hAxes, fac, xLimits, yLimits)

  if (nargin < 3) || isempty(xLimits)
    xLimits = get(hAxes,'XLim');
  end
  if (nargin < 4) || isempty(yLimits)
    yLimits = get(hAxes,'YLim');
  end

  logScale = diff(yLimits)/diff(xLimits);
  powerScale = diff(log10(yLimits))/diff(log10(xLimits));

  set(hAxes,'Xlim',xLimits,...
            'YLim',yLimits,...
            'DataAspectRatio', [1 fac*logScale/powerScale 1]);

end
