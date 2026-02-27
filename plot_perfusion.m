function plot_perfusion(cgrph, u, opts)

arguments
    cgrph chunkgraph
    u     struct
    opts.nplot = 400;
    opts.linewidth = 3;
    opts.contour = false;
end

% Evaluate and plot the solution in the tissue
dom = reshape([min(cgrph) max(cgrph)].', [1 4]);
padx = 0.1*abs(diff(dom(1:2)));
pady = 0.1*abs(diff(dom(3:4)));
dom = dom + [-padx padx -pady pady];
[xx, yy] = meshgrid(linspace(dom(1), dom(2), opts.nplot), ...
                    linspace(dom(3), dom(4), opts.nplot));
eval_opts = [];
eval_opts.accel = true;
eval_opts.forcesmooth = true;
uu = u.bulk(xx, yy, eval_opts);

pcolor(xx, yy, uu, edgecolor='none')
hold on
if ( opts.contour )
    contour(xx, yy, uu, 'k')
end

% Plot the solution in the channel
uc = u.channel;
nedges = length(cgrph.echnks);
xc = cell(nedges, 1);
yc = cell(nedges, 1);
for k = 1:nedges
    xc{k} = cgrph.echnks(k).r(1,:,:); xc{k} = xc{k}(:);
    yc{k} = cgrph.echnks(k).r(2,:,:); yc{k} = yc{k}(:);
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
