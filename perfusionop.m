classdef perfusionop %#ok<*PROP,*PROPLC>

    properties
        cgrph
        root
        xi
        alpha
        phi
        P
        nverts
        nedges
        n
        nch
        npt
        nbc
        ndof
        sigma_starts
        L
        s
        introw1
        introw2
        iintegral
        outer
    end

    methods

        function op = perfusionop(cgrph, opts)
            arguments
                cgrph            chunkgraph
                opts.xi    (1,1) double = 0
                opts.alpha (1,1) double = 0
                opts.phi   (1,1) double = 0
                opts.P     (1,1) double = 0.03
                opts.root  (1,:) double = 1
                opts.outer (1,1) logical = false
            end
            vec = @(x) x(:);
            op.cgrph        = cgrph;
            op.nverts       = length(cgrph.verts);
            op.nedges       = length(cgrph.echnks);
            op.n            = cgrph.k;
            op.nch          = sum([cgrph.echnks.nch]);
            op.npt          = cgrph.npt;
            op.nbc          = 2*op.nedges;
            op.ndof         = op.npt + op.nbc;
            op.sigma_starts = cumsum([1 cgrph.echnks.npt]);
            op.xi           = opts.xi;
            op.alpha        = opts.alpha;
            op.phi          = opts.phi;
            op.P            = opts.P;
            op.root         = opts.root;
            op.L            = arrayfun(@(chnkr) sum(chnkr.wts, 'all'),    cgrph.echnks);
            op.s            = arrayfun(@(chnkr) vec(arclengthfun(chnkr)), cgrph.echnks, uniformoutput=false);
            op.introw1      = arrayfun(@(chnkr) introw(chnkr, 1),         cgrph.echnks, uniformoutput=false);
            op.introw2      = arrayfun(@(chnkr) introw(chnkr, 2),         cgrph.echnks, uniformoutput=false);
            op.iintegral    = @iintegral;

            if ( opts.outer )
                scale = 1.3;
                [center, radius] = minboundcircle(op.cgrph.r(1,:,:), op.cgrph.r(2,:,:));
                cparams = [];
                cparams.nover = 3;
                pref = [];
                pref.k = cgrph.k;
                chnkr = chunkerfunc(@(t) chnk.curves.bymode(t, scale*radius, center), cparams, pref);
                op.outer = chnkr;
                op.ndof = op.ndof + op.outer.npt;
            end

        end

        function varargout = size(op, dim)
            [m, n] = deal(op.ndof);
            if ( nargout < 2 )
                if ( nargin == 1 )
                    varargout = {[m, n]};
                elseif ( dim == 1 )
                    varargout = {m};
                elseif ( dim == 2 )
                    varargout = {n};
                end
            else
                varargout = {m, n};
            end

        end

        function [A, rhs, R] = discretize(op, opts)
            arguments
                op
                opts.blocks  = false
                opts.fmm     = false
                opts.precond = false
            end

            if ( opts.fmm )
                opts.blocks = true;
            end

            if ( isempty(op.outer) )
                A = cell(2);
                R = cell(2);
            else
                A = cell(3);
                R = cell(3);
            end

            [xi, alpha, phi, P] = deal(op.xi, op.alpha, op.phi, op.P);

            Ls = op.L;
            introw1 = op.introw1;
            introw2 = op.introw2;

            if ( opts.fmm )
                cgrph = op.cgrph;
                chnkrs = merge(op.cgrph.echnks);
                Skern = kernel('helm', 's', 1i*phi);
                cormat = chunkermat(op.cgrph, Skern, struct(rcip=false, corrections=true));
                S = @(sigma,tol) real(chunkermatapply(chnkrs, Skern, sigma, cormat, ...
                    struct(rcip=false, accel=true, eps=tol)));
                A{1,1} = @(sigma,tol) S(sigma,tol) + alpha*sigma - P*xi*iintegral(cgrph, sigma, 2);

                if ( opts.precond )
                    Dpkern = kernel('helm', 'dprime', 1i*phi);
                    %cormat = chunkermat(op.cgrph, Dpkern, struct(rcip=false, corrections=true));
                    %Dp = @(sigma,tol) real(chunkermatapply(chnkrs, Dpkern, sigma, cormat, ...
                    %    struct(rcip=false, accel=true, eps=tol)));
                    %R{1,1} = @(sigma,tol) Dp(sigma,tol);
                    chopts = [];
                    chopts.rcip = false;
                    Dps = arrayfun(@(c) chunkermat(c, Dpkern, chopts), op.cgrph.echnks, uniformoutput=false);
                    Dp = matlab.internal.math.blkdiag(Dps{:});
                    R{1,1} = @(sigma,tol) Dp*sigma;
                else
                    R{1,1} = @(sigma,tol) sigma;
                end
            else
                chopts = [];
                chopts.rcip = false;
                Skern = kernel('helm', 's', 1i*phi);
                S = chunkermat(op.cgrph, Skern, chopts);
                I = eye(size(S));
                Qs  = arrayfun(@myintmat, op.cgrph.echnks, uniformoutput=false);
                Qs2 = cellfun(@(x) x^2, Qs, uniformoutput=false);
                Q2 = matlab.internal.math.blkdiag(Qs2{:});
                A{1,1} = S + alpha*I - P*xi*Q2;
                if ( opts.precond )
                    Dpkern = kernel('helm', 'dprime', 1i*phi);
                    %Dp = chunkermat(op.cgrph, Dpkern, chopts);
                    Dps  = arrayfun(@(c) chunkermat(c, Dpkern, chopts), op.cgrph.echnks, uniformoutput=false);
                    Dp = matlab.internal.math.blkdiag(Dps{:});
                    R{1,1} = Dp;
                else
                    R{1,1} = speye(op.npt);
                end
            end

            s = matlab.internal.math.blkdiag(op.s{:});
            one = spones(s);
            A{1,2} = [-P*one  -P*s];

            nedges = op.nedges;
            npt = op.npt;
            isroot = false(op.nverts, 1);
            isroot(op.root) = true;
            bcval = zeros(op.nbc, 1);

            degs = cellfun(@(x) length(x{1}), op.cgrph.vstruc);
            assert(sum(degs) == op.nbc);
            bc = zeros(op.nbc, op.npt+op.nbc);
            ibc = 1;
            sigma_starts = op.sigma_starts;

            % Note: the following does not support self loops, but could be modified to
            for k_vert = 1:op.nverts
                edge_ids = op.cgrph.vstruc{k_vert}{1};
                edge_dir = op.cgrph.vstruc{k_vert}{2};
                deg = degs(k_vert);
                % Impose deg-1 pairwise continuity conditions
                for k_edge = 1:length(edge_ids)-1
                    edge_id1 = edge_ids(k_edge);
                    edge_id2 = edge_ids(k_edge+1);
                    sigma1 = zeros(1, npt);
                    sigma2 = zeros(1, npt);
                    a1 = zeros(1, nedges);
                    a2 = zeros(1, nedges);
                    b1 = zeros(1, nedges);
                    b2 = zeros(1, nedges);
                    a1(edge_id1) = 1;
                    a2(edge_id2) = 1;
                    if edge_dir(k_edge) == 1
                        idx = sigma_starts(edge_id1):sigma_starts(edge_id1+1)-1;
                        sigma1(idx) = xi * introw2{edge_id1};
                        b1(edge_id1) = Ls(edge_id1);
                    end
                    if edge_dir(k_edge+1) == 1
                        idx = sigma_starts(edge_id2):sigma_starts(edge_id2+1)-1;
                        sigma2(idx) = xi * introw2{edge_id2};
                        b2(edge_id2) = Ls(edge_id2);
                    end
                    sigma_cont = sigma1 - sigma2;
                    a_cont     = a1 - a2;
                    b_cont     = b1 - b2;
                    bc(ibc,:) = [sigma_cont a_cont b_cont];
                    ibc = ibc+1;
                end

                % Impose 1 Kirchoff condition
                sigma_kirchoff = zeros(1, op.cgrph.npt);
                a_kirchoff     = zeros(1, nedges);
                b_kirchoff     = zeros(1, nedges);
                if (isroot(k_vert))
                    k_edge = 1;
                    edge_id = edge_ids(k_edge);
                    a_kirchoff(edge_id) = 1;
                    if edge_dir(k_edge) == 1
                        idx = sigma_starts(edge_id):sigma_starts(edge_id+1)-1;
                        sigma_kirchoff(idx) = xi * introw2{edge_id};
                        b_kirchoff(edge_id) = Ls(edge_id);
                    end
                    bcval(ibc) = 1;
                else
                    for k_edge = 1:length(edge_ids)
                        edge_id = edge_ids(k_edge);
                        idx = sigma_starts(edge_id):sigma_starts(edge_id+1)-1;
                        % Note: The flux on each edge is given by 1/xi * U'.
                        % But here, we scale the flux by an overall factor of xi
                        % to gracefully handle the case when xi = 0.
                        if edge_dir(k_edge) == 1
                            sigma_kirchoff(idx) = xi * introw1{edge_id};
                            b_kirchoff(edge_id) = 1;
                        elseif edge_dir(k_edge) == -1
                            b_kirchoff(edge_id) = -1;
                        end
                    end
                end
                bc(ibc,:) = [sigma_kirchoff a_kirchoff b_kirchoff];
                ibc = ibc+1;
            end
            bc = sparse(bc);

            A{2,1} = bc(:,1:op.npt);
            A{2,2} = bc(:,op.npt+1:end);

            R{2,2} = speye(op.nbc);
            R{1,2} = sparse(op.npt, op.nbc);
            R{2,1} = sparse(op.nbc, op.npt);

            rhs = cell(2, 1);
            rhs{1} = zeros(op.npt, 1);
            rhs{2} = bcval;

            if ( ~isempty(op.outer) )

                if ( opts.fmm )
                    outer = op.outer;
                    echnks = merge(op.cgrph.echnks);
                    Skern  = kernel('helm', 's',      1i*phi);
                    Spkern = kernel('helm', 'sprime', 1i*phi);

                    cormat = chunkermat(op.outer, Spkern, struct(rcip=false, corrections=true));
                    Sp = @(tau,tol) real(chunkermatapply(outer, Spkern, tau, cormat, ...
                        struct(rcip=false, accel=true, eps=tol)));
                    A{3,3} = @(tau,tol) 0.5*tau + Sp(tau,tol);
                    R{3,3} = @(tau,tol) tau;

                    cormat = chunkerkernevalmat(op.outer, Skern, op.cgrph, struct(corrections=true));
                    out2in = @(tau,tol) real(chunkerkerneval(outer, Skern, tau, echnks, ...
                        struct(cormat=cormat, accel=true, eps=tol)));

                    cormat = chunkerkernevalmat(op.cgrph, Spkern, op.outer, struct(corrections=true));
                    in2out = @(sigma,tol) real(chunkerkerneval(echnks, Spkern, sigma, outer, ...
                        struct(cormat=cormat, accel=true, eps=tol)));
                else
                    chopts = [];
                    chopts.rcip = false;
                    Skern  = kernel('helm', 's',      1i*phi);
                    Spkern = kernel('helm', 'sprime', 1i*phi);
                    Sp = chunkermat(op.outer, Spkern, chopts);
                    I = eye(size(Sp));
                    A{3,3} = 0.5*I + Sp;
                    R{3,3} = speye(op.outer.npt);
                    out2in = chunkerkernevalmat(op.outer, Skern,  op.cgrph);
                    in2out = chunkerkernevalmat(op.cgrph, Spkern, op.outer);
                end

                A{1,3} = out2in;
                A{3,1} = in2out;
                A{2,3} = sparse(op.nbc, op.outer.npt);
                A{3,2} = sparse(op.outer.npt, op.nbc);

                R{1,3} = sparse(op.npt, op.outer.npt);
                R{2,3} = sparse(op.nbc, op.outer.npt);
                R{3,1} = sparse(op.outer.npt, op.npt);
                R{3,2} = sparse(op.outer.npt, op.nbc);

                rhs{3} = zeros(op.outer.npt, 1);
            end

            if ( ~opts.blocks )
                % Collapse output into a dense matrix
                A = cell2mat(cellfun(@full, A, uniformoutput=false));
                R = cell2mat(cellfun(@full, R, uniformoutput=false));
                rhs = cell2mat(rhs);
            end

        end

        function [u, vars] = solve(op, opts)

            arguments
                op
                opts.method  = 'auto'
                opts.fmm     = false
                opts.schur   = false
                opts.precond = false
                opts.tol     = 1e-10
                opts.maxit   = op.ndof
            end

            if ( strcmpi(opts.method, 'auto') )
                if ( op.ndof < 10000 )
                    opts.method = 'direct';
                    opts.fmm    = false;
                    opts.schur  = false;
                else
                    opts.method = 'gmres';
                    opts.fmm    = true;
                    opts.schur  = true;
                end
            end

            if ( strcmpi(opts.method, 'direct') && opts.fmm )
                warning('Direct solve cannot use FMM.');
                opts.fmm = false;
            end

            switch lower(opts.method)
                case 'direct'
                    solve = @(A,b,~) A \ b;
                case 'gmres'
                    solve = @(A,b,R) right_gmres(A, b, [], opts.tol, min(opts.maxit, length(b)), R);
                otherwise
                    error('Unknown method ''%s''', opts.method);
            end

            [A, rhs, R] = discretize(op, fmm=opts.fmm, blocks=opts.schur, precond=opts.precond);

            if ( ~opts.schur )
                if ( opts.fmm )
                    ii = 1:op.npt;
                    bb = op.npt+(1:op.nbc);
                    tol = opts.tol;
                    tol = 1e-10;
                    if ( ~isempty(op.outer) )
                        jj = (op.npt+op.nbc+1):op.ndof;
                        A = @(x) [ A{1,1}(x(ii),tol) + A{1,2}*x(bb) + A{1,3}(x(jj),tol)  ;
                                   A{2,1}*x(ii)      + A{2,2}*x(bb)                      ;
                                   A{3,1}(x(ii),tol)                + A{3,3}(x(jj),tol) ];
                        R = @(x) [ R{1,1}(x(ii),tol) ; R{2,2}*x(bb) ; R{3,3}*x(jj) ];
                        rhs = [rhs{1} ; rhs{2} ; rhs{3}];
                    else
                        A = @(x) [ A{1,1}(x(ii),tol) + A{1,2}*x(bb) ;
                                   A{2,1}*x(ii)      + A{2,2}*x(bb) ];
                        R = @(x) [ R{1,1}(x(ii),tol) ; R{2,2}*x(bb) ];
                        rhs = [rhs{1} ; rhs{2}];
                    end
                end
                sigma_a_b = solve(A, rhs, R);
            else
                Abb = decomposition(A{2,2});
                rhsb = rhs{2};
                if ( ~isempty(op.outer) )
                    if ( opts.fmm )
                        ii = 1:op.npt;
                        jj = op.npt+(1:op.outer.npt);
                        Aii = @(x,tol) [ A{1,1}(x(ii),tol) + A{1,3}(x(jj),tol) ;
                                         A{3,1}(x(ii),tol) + A{3,3}(x(jj),tol) ];
                        Rii = @(x,tol) [ R{1,1}(x(ii),tol) ; R{3,3}(x(jj),tol) ];
                    else
                        Aii = [ A{1,1} A{1,3} ;
                                A{3,1} A{3,3} ];
                        Rii = [ R{1,1} R{1,3} ;
                                R{3,1} R{3,3} ];
                    end
                    Aib = [A{1,2}; A{3,2}];
                    Abi = [A{2,1}  A{2,3}];
                    rhsi = [rhs{1} ; rhs{3}];
                else
                    Aii = A{1,1};
                    Aib = A{1,2};
                    Abi = A{2,1};
                    Rii = R{1,1};
                    rhsi = rhs{1};
                end
                if ( opts.fmm )
                    tol = opts.tol;
                    tol = 1e-10;
                    Sii = @(x) Aii(x,tol) - Aib*(Abb\(Abi*x));
                    Rii = @(x) Rii(x,tol);
                else
                    Sii = Aii - Aib*(Abb\Abi);
                end
                fi = rhsi - Aib*(Abb\rhsb);
                sigma = solve(Sii, fi, Rii);
                a_b = Abb \ (rhsb-Abi*sigma);
                if ( ~isempty(op.outer) )
                    tau   = sigma(op.npt+1:end);
                    sigma = sigma(1:op.npt);
                    sigma_a_b = [sigma ; a_b ; tau];
                else
                    sigma_a_b = [sigma ; a_b];
                end
            end

            sigma = sigma_a_b(1:op.npt);
            a     = sigma_a_b(op.npt+1:op.npt+op.nedges);
            b     = sigma_a_b(op.npt+op.nedges+1:op.npt+2*op.nedges);
            if ( ~isempty(op.outer) )
                tau = sigma_a_b(op.npt+2*op.nedges+1:end);
            end

            u = [];
            u.channel = cell(op.nedges, 1);
            for k = 1:op.nedges
                idx = op.sigma_starts(k):op.sigma_starts(k+1)-1;
                sigma_k = sigma(idx);
                aa = repmat(a(k), length(sigma_k), 1);
                bb = repmat(b(k), length(sigma_k), 1);
                chnkr_k = op.cgrph.echnks(k);
                u.channel{k} = op.xi*iintegral(chnkr_k, sigma_k, 2) + aa + bb.*op.s{k};
            end

            cgrph = op.cgrph;
            outer = op.outer;
            Skern = kernel('helm', 's', 1i*op.phi);
            function u = u_bulk(x, y, opts)
                if ( nargin < 3 )
                    opts = [];
                end
                u = chunkerkerneval(cgrph, Skern, sigma, [x(:) y(:)].', opts);
                if ( ~isempty(outer) )
                    u = u + chunkerkerneval(outer, Skern, tau, [x(:) y(:)].', opts);
                end
                u = reshape(u, size(x));
                u = real(u);
            end
            u.bulk = @u_bulk;

            vars = [];
            vars.sigma = sigma;
            vars.a = a;
            vars.b = b;
            if ( ~isempty(op.outer) )
                vars.tau = tau;
            end

        end

        function check_conservation(op, u, vars, opts)

            arguments
                op
                u    struct
                vars struct
                opts.ngrid (1,1) double = 800
            end

            fprintf('Values:\n\n');

            % idx = op.sigma_starts(iedge):op.sigma_starts(iedge+1)-1;
            % Qsigma = op.iintegral(op.cgrph.echnks(iedge), sigma(idx));
            % Qch = Qsigma((ich-1)*op.n + (1:op.n));
            % flux_root = -1/op.xi * (interp_start * Qch + b(iedge));

            flux_root1 = -1/op.xi * vars.b(op.root);
            fprintf('    (1)  (Qσ(root)+b)/ξ  = %.16g\n', flux_root1);

            iedge = 1;
            ich = 1;
            Dleg = lege.dermat(op.n);
            ds = vecnorm(op.cgrph.echnks(iedge).d(:,:,ich)).';
            Dch = Dleg ./ ds;
            interp_start = lege.matrin(op.n, -1);
            interp_end   = lege.matrin(op.n,  1);
            flux_root2 = -1/op.xi * interp_start * Dch * u.channel{iedge}((ich-1)*op.n + (1:op.n));
            fprintf('    (2)  U''(root)/ξ      = %.16g\n', flux_root2);

            w = op.cgrph.wts(:);
            int_sigma = sum(vars.sigma .* w);
            fprintf('    (3)  ∫σ              = %.16g\n', int_sigma);

            if ( ~isempty(op.outer) )
                dom = reshape([min(op.outer) max(op.outer)].', [1 4]);
                [xx, yy] = meshgrid(linspace(dom(1), dom(2), opts.ngrid), ...
                                    linspace(dom(3), dom(4), opts.ngrid));
                in = chunkerinterior(op.outer, [xx(:) yy(:)].');
                in = reshape(in, size(xx));
            else
                % As r -> inf, K0(phi*r) -> sqrt(pi/(2*phi*r))*exp(-phi*r).
                % So K0 ~ 1e-16 when r ~ 35/phi.
                dom = reshape([min(op.cgrph) max(op.cgrph)].', [1 4]);
                pad = 35/op.phi;
                dom = dom + [-pad pad -pad pad];
                [xx, yy] = meshgrid(linspace(dom(1), dom(2), opts.ngrid), ...
                                    linspace(dom(3), dom(4), opts.ngrid));
                in = true(size(xx));
            end
            eval_opts = [];
            eval_opts.accel = true;
            eval_opts.tol = 1e-10;
            eval_opts.forcesmooth = true;
            uu = nan(size(xx));
            uu(in) = u.bulk(xx(in), yy(in), eval_opts);

            dx = diff(dom(1:2))/(opts.ngrid-1);
            dy = diff(dom(3:4))/(opts.ngrid-1);
            int_u = op.phi^2 * sum(uu(in), 'all') * dx * dy;
            fprintf('    (4)  ∫uϕ²            = %.16g\n', int_u);

            fprintf('\nConservation relative error:\n\n');
            fprintf('    |(1)-(4)| / |(1)| = %g\n', abs(flux_root1 - int_u) / abs(flux_root1));
            fprintf('\n');
        end

    end

end

function imat = myintmat(chnkr)
%INTMAT returns the matrix of integration in arc length along the 
% chunker object. The constant of integration is selected such that the
% integral is zero at the left end of the first chunk of the chunker. 

    n = chnkr.k;
    nch = chnkr.nch;
	panel_int = lege.intmat(n);
	ds = reshape(vecnorm(chnkr.d), [1 n nch]);
	panel_int = panel_int .* ds;
	[~, w, ~, ~] = lege.exps(n);
	ds = reshape(w, [n 1]) .* reshape(ds, [n nch]);
	% Set diagonal blocks
	panel_int = num2cell(panel_int, [1 2]);
	imat = blkdiag(panel_int{:});
	% maybe not the fast way to do this, but find the panel ordering
	order_vec = zeros(nch, 1);
	j = 1;
	for i = 1:nch
		order_vec(i) = j;
		j = chnkr.adj(2, j);
	end
	% Set off diagonal matrix components
    for i = 2:nch
		for j = 1:i-1
			u = 1 + (order_vec(i)-1)*n;
			v = 1 + (order_vec(j)-1)*n;
			tmp = ones(n, 1) * ds(:, order_vec(j)).';
			imat(u:u+n-1, v:v+n-1) = tmp;
		end
    end

end

function int = dintegral(chnkr, f)
    int = sum(chnkr.wts(:) .* f);
end

function int = iintegral(chnkr, f, p)

    if ( nargin < 3 )
        p = 1;
    end

    if ( isa(chnkr, 'chunkgraph') )
        cgrph = chnkr;
        int = zeros(cgrph.npt, 1);
        starts = cumsum([1 cgrph.echnks.npt]);
        for k = 1:length(cgrph.echnks)
            chnkr = cgrph.echnks(k);
            idx = starts(k):starts(k+1) - 1;
            int(idx) = iintegral(chnkr, f(idx), p);
        end
        return
    end

    n   = chnkr.k;
    nch = chnkr.nch;

    persistent Q
    if ( isempty(Q) || size(Q,1) ~= n )
        Q = lege.intmat(n);
    end

    ds = reshape(arclengthdens(chnkr), [1 n nch]);
    imat = Q .* ds;
    w = reshape(chnkr.wts, [n 1 nch]);
    int = reshape(f, [n 1 nch]);
    for j = 1:p
        sums = sum(w.*int);
        sums = sums(:,:,1:end-1);
        int = pagemtimes(imat, int);
        int = int + cumsum(cat(3, 0, sums));
    end
    int = reshape(int, size(f));

end

function int = iintegralT(chnkr, f, p)

    if ( nargin < 3 )
        p = 1;
    end

    if ( isa(chnkr, 'chunkgraph') )
        cgrph = chnkr;
        int = zeros(cgrph.npt, 1);
        starts = cumsum([1 cgrph.echnks.npt]);
        for k = 1:length(cgrph.echnks)
            chnkr = cgrph.echnks(k);
            idx = starts(k):starts(k+1) - 1;
            int(idx) = iintegralT(chnkr, f(idx), p);
        end
        return
    end

    n   = chnkr.k;
    nch = chnkr.nch;

    persistent Q
    if ( isempty(Q) || size(Q,1) ~= n )
        Q = lege.intmat(n);
    end

    ds = reshape(arclengthdens(chnkr), [1 n nch]);
    imatT = pagetranspose(Q .* ds);
    w = reshape(chnkr.wts, [n 1 nch]);
    int = reshape(f, [n 1 nch]);
    for j = 1:p
        sums = flip(cumsum(flip(sum(int))));
        sums = sums(:,:,2:end);
        int = pagemtimes(imatT, int);
        int = int + w.*cat(3, sums, 0);
    end
    int = reshape(int, size(f));

end

function introw = introw(chnkr, p)

    if ( nargin < 2 )
        p = 1;
    end

    w = chnkr.wts(:);
    introw = iintegralT(chnkr, w, p-1).';

end

function [x, flag, relres, iter, resvec] = right_gmres(A, b, restart, tol, maxit, P)
%RIGHT_GMRES   GMRES with right preconditioning.

if ( isa(A, 'function_handle') )
    afun = A;
else
    afun = @(x) A*x;
end

if ( isa(P, 'function_handle') )
    pfun = P;
else
    pfun = @(x) P*x;
end

% Define new linear operator A*P
AP = @(x) afun(pfun(x));

% Solve system A*P*y = b
[y, flag, relres, iter, resvec] = gmres(AP, b, restart, tol, maxit);

fprintf('gmres converged at iteration %d to a solution with relative residual %g.\n', iter(2), relres);

% Define x = P*y
x = pfun(y);

end

function B = blockdiag(A, n)

m = size(A, 1);
p = m / n;
[ii, jj] = ndgrid(1:n);
kk = reshape(1:p, [1 1 p]);
i = reshape(ii + n*(kk-1), [], 1);
j = reshape(jj + n*(kk-1), [], 1);
idx = sub2ind([m m], i, j);
v = A(idx);
B = sparse(i, j, v, m, m);

end
