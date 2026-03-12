function plot_perfusion(op, u, opts)

arguments
    op perfusionop
    u  struct
    opts.nplot = 400;
    opts.linewidth = 3;
    opts.contour = false;
    opts.partialpressure = true;
    opts.forcesmooth = true;
end

% Evaluate and plot the solution in the tissue
if ( ~isempty(op.outer) )
    dom = reshape([min(op.outer) max(op.outer)].', [1 4]);
    [xx, yy] = meshgrid(linspace(dom(1), dom(2), opts.nplot), ...
                        linspace(dom(3), dom(4), opts.nplot));
    in = chunkerinterior(op.outer, [xx(:) yy(:)].');
    in = reshape(in, size(xx));
else
    dom = reshape([min(op.cgrph) max(op.cgrph)].', [1 4]);
    padx = 0.1*abs(diff(dom(1:2)));
    pady = 0.1*abs(diff(dom(3:4)));
    pad = max(padx, pady);
    padx = max(pad, padx);
    pady = max(pad, pady);
    dom = dom + [-padx padx -pady pady];
    [xx, yy] = meshgrid(linspace(dom(1), dom(2), opts.nplot), ...
                        linspace(dom(3), dom(4), opts.nplot));
    in = true(size(xx));
end
eval_opts = [];
eval_opts.accel = true;
eval_opts.tol = 1e-10;
eval_opts.forcesmooth = opts.forcesmooth;
uu = nan(size(xx));
uu(in) = u.bulk(xx(in), yy(in), eval_opts);

if ( opts.partialpressure )
    title('Partial pressure')
    uu = uu/op.P;
else
    title('Concentration')
end
hold on
pcolor(xx, yy, uu, edgecolor='none')
shading interp
if ( opts.contour )
    contour(xx, yy, uu, 'k')
end

% Plot the solution in the channel
uc = u.channel;
nedges = length(op.cgrph.echnks);
xc = cell(nedges, 1);
yc = cell(nedges, 1);
for k = 1:nedges
    xc{k} = op.cgrph.echnks(k).r(1,:,:); xc{k} = xc{k}(:);
    yc{k} = op.cgrph.echnks(k).r(2,:,:); yc{k} = yc{k}(:);
end
nans = repmat({NaN}, nedges, 1);
xc = [xc nans].'; xc = cat(1, xc{:});
yc = [yc nans].'; yc = cat(1, yc{:});
uc = [uc nans].'; uc = cat(1, uc{:});

patch(xc, yc, uc, edgecolor='interp', linewidth=opts.linewidth)
hold off

axis equal
axis(dom)
clim([0 1])
colorbar

end
