% Compare GMRES convergence with and without right preconditioning across
% inverse permeabilities lambda.

%% Graph parameters
name = 'L2_1_Tr9L';   % Terminal cell (see data/)
smooth = 1;
fittol = 1e-4;
n = 16;
nover = 1;
cgrph = graphtools.fit(graphtools.gallery(name), smooth=smooth, tol=fittol, order=n, nover=nover);

%% Physical parameters
root   = 1;     % Indices of trachea roots
alpha  = 250;   % Scaled channel area (1/xi)
phi    = 0.1;   % Yukawa parameter
eta    = 0.03;  % Henry's constant
outer  = true;

%%% Solver parameters
method   = 'gmres';
schur    = true;
fmm      = false;
tol      = 1e-10;
lambda = [0 0.01 0.1 1];

resvec_pre   = cell(length(lambda), 1);
resvec_unpre = cell(length(lambda), 1);
for k = 1:length(lambda)
    fprintf('lambda = %g, precond = true\n', lambda(k));
    op = perfusionop(cgrph, root=root, outer=outer, alpha=alpha, lambda=lambda(k), phi=phi, eta=eta);
    [~, ~, info] = op.solve(method=method, schur=schur, fmm=fmm, precond=true, tol=tol);
    resvec_pre{k} = info.resvec / info.resvec(1);
    fprintf('lambda = %g, precond = false\n', lambda(k));
    [~, ~, info] = op.solve(method=method, schur=schur, fmm=fmm, precond=false, tol=tol);
    resvec_unpre{k} = info.resvec / info.resvec(1);
end

%% Plot
% figure(1), clf
% for k = 1:length(lambda)
%     semilogy(resvec_pre{k}, '-o')
%     hold on
% end
% hold off
% grid on
% title('Preconditioned')
% legend('\lambda = 0', '\lambda = 0.01', '\lambda = 0.1', '\lambda = 1')
% 
% figure(2), clf
% for k = 1:length(lambda)
%     semilogy(resvec_unpre{k}, '-o')
%     hold on
% end
% hold off
% grid on
% title('Unpreconditioned')
% legend('\lambda = 0', '\lambda = 0.01', '\lambda = 0.1', '\lambda = 1')

close all
figure(1)
lw = 2;
%cmap = redblue(5);
values = [0.8 0.82 0.99; 0 0 1];
m = 5;
P = size(values, 1);
cmap = interp1(1:P, values, linspace(1,P,m), 'linear');
cmap = colororder;
%cmap = nebula(50);
%cmap = cmap([1 13 20 45],:);
cmap = nebula(4);
cmap = flip(cmap);

semilogy(resvec_unpre{4}, '-',  linewidth=lw, color=cmap(4,:))
hold on
semilogy(resvec_unpre{3}, '-',  linewidth=lw, color=cmap(3,:))
semilogy(resvec_pre{2},   ':', linewidth=lw, color=cmap(2,:))
semilogy(resvec_pre{1},   ':', linewidth=lw, color=cmap(1,:))
semilogy(resvec_unpre{1}, '-',  linewidth=lw, color=cmap(1,:))

hold off
grid on
ylim([1e-10 1])
axis square
xlabel('GMRES iteration')
ylabel('Relative residual')
%legend('\lambda = 1', '\lambda = 0.1', '\lambda = 0.01 (precond)',  '\lambda = 0 (precond)', '\lambda = 0 (no precond)')
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.FontSize = 12;
ax.LineWidth = 1;
ax.GridLineWidth = 0.7;

colormap(nebula(300))
cb = colorbar(location='southoutside');
%ticks = ([1 13 20 45]-1)/49;
ticks = linspace(0, 1, 4);
cb.Limits = [0 1];
cb.Ticks = ticks;
cb.TickLabels = flip({'$\lambda = 0$', '$\lambda = 10^{-2}$', '$\lambda = 10^{-1}$', '$\lambda = 1$'});
cb.TickLabelInterpreter = 'latex';
cb.FontSize = 12;
cb.AxisLocation = 'in';
cb.LineWidth = 1;
pos = cb.Position;
pos(1) = pos(1) + 0.22;
pos(2) = pos(2) + 0.72;
pos(3) = pos(3) - 0.15;
pos(4) = pos(4) - 0.005;
cb.Position = pos;
cb.TickDirection = 'out';
cb.TickLength = 0.02;
cb.Ruler.TickLabelRotation = 0;

shg
