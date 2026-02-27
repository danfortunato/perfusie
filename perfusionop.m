classdef perfusionop %#ok<*PROP,*PROPLC>

    properties
        cgrph
        root
        xi
        alpha
        phi

        nverts
        nedges
        n
        nch
        npt
        nbc
        ndof
        sigma_starts

        s
        L
        introw1
        introw2
    end

    methods

        function op = perfusionop(cgrph, opts)
            arguments
                cgrph            chunkgraph
                opts.xi    (1,1) double = 0
                opts.alpha (1,1) double = 0
                opts.phi   (1,1) double = 0
                opts.root  (1,:) double = 1
            end
            op.cgrph  = cgrph;
            op.nverts = length(cgrph.verts);
            op.nedges = length(cgrph.echnks);
            op.n      = cgrph.k;
            op.nch    = sum([cgrph.echnks.nch]);
            op.npt    = cgrph.npt;
            op.nbc    = 2*op.nedges;
            op.ndof   = op.npt + op.nbc;
            op.sigma_starts = cumsum([1 cgrph.echnks.npt]);

            op.xi    = opts.xi;
            op.alpha = opts.alpha;
            op.phi   = opts.phi;
            op.root  = opts.root;

            op.L = arrayfun(@(chnkr) sum(chnkr.wts, 'all'), cgrph.echnks);
            op.s = arclengthfun_graph(cgrph);
            op.introw1 = arrayfun(@(chnkr) introw(chnkr, 1), cgrph.echnks, uniformoutput=false);
            op.introw2 = arrayfun(@(chnkr) introw(chnkr, 2), cgrph.echnks, uniformoutput=false);

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

        function [A, rhs] = discretize(op, opts)
            arguments
                op
                opts.blocks = false
                opts.fmm  = false
            end

            if ( opts.fmm )
                opts.blocks = true;
            end

            A = cell(2);

            [xi, alpha, phi] = deal(op.xi, op.alpha, op.phi);

            Ls = op.L;
            %Qs  = op.Q;
            %Qs2 = op.Q2;
            %interp_end = op.interp_end;
            introw1 = op.introw1;
            introw2 = op.introw2;

            if ( opts.fmm )
                cgrph = op.cgrph;
                chnkrs = merge(op.cgrph.echnks);
                Skern = kernel('helm', 's', 1i*phi);
                cormat = chunkermat(op.cgrph, Skern, struct(rcip=false, corrections=true));
                S = @(sigma,tol) real(chunkermatapply(chnkrs, Skern, sigma, cormat, ...
                    struct(rcip=false, accel=true, eps=tol)));
                %A{1,1} = @(sigma,tol) S(sigma,tol) + alpha*sigma - xi*Q2*sigma;
                A{1,1} = @(sigma,tol) S(sigma,tol) + alpha*sigma - xi*iintegral(cgrph, sigma, 2);
            else
                chopts = [];
                chopts.rcip = false;
                chopts.nonsmoothonly = false;
                Skern = kernel('helm', 's', 1i*phi);
                S = chunkermat(op.cgrph, Skern, chopts);
                I = eye(size(S));
                Qs  = arrayfun(@myintmat, op.cgrph.echnks, uniformoutput=false);
                Qs2 = cellfun(@(x) x^2, Qs, uniformoutput=false);
                Q2 = matlab.internal.math.blkdiag(Qs2{:});
                A{1,1} = S + alpha*I - xi*Q2;
            end

            one = edgewise(op.cgrph, @(chnkr) ones(chnkr.npt, 1));
            s   = edgewise(op.cgrph, @(chnkr) reshape(arclengthfun(chnkr), chnkr.npt, 1));
            A{1,2} = [-one  -s];

            nedges = op.nedges;
            npt = op.npt;
            isroot = false(op.nverts, 1);
            isroot(op.root) = true;
            bcval = zeros(op.nbc, 1);

            degs = cellfun(@(x) length(x{1}), op.cgrph.vstruc);
            assert(sum(degs) == op.nbc);
            bc = zeros(op.nbc, op.ndof);
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
                        %b1(edge_id1) = 1;
                    end
                    if edge_dir(k_edge+1) == 1
                        idx = sigma_starts(edge_id2):sigma_starts(edge_id2+1)-1;
                        sigma2(idx) = xi * introw2{edge_id2};
                        b2(edge_id2) = Ls(edge_id2);
                        %b2(edge_id2) = 1;
                    end
                    % sigma1 = sigma1 ./ Ls(edge_id1); sigma2 = sigma2 ./ Ls(edge_id2);
                    % a1 = a1 ./ Ls(edge_id1);         a2 = a2 ./ Ls(edge_id2);
                    % b1 = b1 ./ Ls(edge_id1);         b2 = b2 ./ Ls(edge_id2);
                    % maxL = max(Ls([edge_id1 edge_id2]));
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
                        %b_kirchoff(edge_id) = 1;
                    end
                    bcval(ibc) = 1;
                else
                    for k_edge = 1:length(edge_ids)
                        edge_id = edge_ids(k_edge);
                        idx = sigma_starts(edge_id):sigma_starts(edge_id+1)-1;
                        if edge_dir(k_edge) == 1
                            sigma_kirchoff(idx) = xi * introw1{edge_id};
                            b_kirchoff(edge_id) = 1;
                            %b_kirchoff(edge_id) = 1 / Ls(edge_id);
                        elseif edge_dir(k_edge) == -1
                            b_kirchoff(edge_id) = -1;
                            %b_kirchoff(edge_id) = -1 / Ls(edge_id);
                        end
                    end
                end
                bc(ibc,:) = [sigma_kirchoff a_kirchoff b_kirchoff];
                ibc = ibc+1;
            end
            bc = sparse(bc);

            A{2,1} = bc(:,1:op.npt);
            A{2,2} = bc(:,op.npt+1:end);

            rhs = cell(2, 1);
            rhs{1} = zeros(op.npt, 1);
            rhs{2} = bcval;

            if ( ~opts.blocks )
                % Collapse output into a dense matrix
                A   = [ A{1,1}       full(A{1,2}) ;
                        full(A{2,1}) full(A{2,2}) ];
                rhs = [rhs{1} ; rhs{2}];
            end

        end

        function u = solve(op, opts)

            arguments
                op
                opts.method = 'auto'
                opts.fmm = false
                opts.schur = false
                opts.tol = 1e-10
                opts.maxit = op.ndof
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

            if ( strcmpi(opts.method, 'direct') )
                opts.fmm = false;
            end

            switch lower(opts.method)
                case 'direct'
                    solve = @(A,b) A \ b;
                case 'gmres'
                    solve = @(A,b) gmres(A, b, [], opts.tol, min(opts.maxit, length(b)));
            end

            [A, rhs] = discretize(op, fmm=opts.fmm, blocks=opts.schur);

            if ( ~opts.schur )
                if ( opts.fmm )
                    ii = 1:op.npt;
                    bb = op.npt+(1:op.nbc);
                    tol = opts.tol;
                    A = @(x) [ A{1,1}(x(ii),tol) + A{1,2}*x(bb) ;
                               A{2,1}*x(ii)      + A{2,2}*x(bb) ];
                    rhs = [rhs{1} ; rhs{2}];
                end
                sigma_a_b = solve(A, rhs);
            else
                Abb = decomposition(A{2,2});
                if ( opts.fmm )
                    tol = opts.tol;
                    Sii = @(x) A{1,1}(x,tol) - A{1,2}*(Abb\(A{2,1}*x));
                else
                    Sii = A{1,1} - A{1,2}*(Abb\A{2,1});
                end
                fi = rhs{1} - A{1,2}*(Abb\rhs{2});
                sigma = solve(Sii, fi);
                a_b = Abb \ (rhs{2}-A{2,1}*sigma);
                sigma_a_b = [sigma ; a_b];
            end

            sigma = sigma_a_b(1:op.npt);
            a     = sigma_a_b(end-2*op.nedges+1:end-op.nedges);
            b     = sigma_a_b(end-op.nedges+1:end);

            u = [];
            u.channel = cell(op.nedges, 1);
            for k = 1:op.nedges
                idx = op.sigma_starts(k):op.sigma_starts(k+1)-1;
                sigma_k = sigma_a_b(idx);
                aa = repmat(a(k), length(sigma_k), 1);
                bb = repmat(b(k), length(sigma_k), 1);
                chnkr_k = op.cgrph.echnks(k);
                u.channel{k} = op.xi*iintegral(chnkr_k, sigma_k, 2) + aa + bb.*op.s{k};
            end

            cgrph = op.cgrph;
            kern = kernel('helm', 's', 1i*op.phi);
            function u = u_bulk(x, y, opts)
                if ( nargin < 3 )
                    opts = [];
                end
                u = chunkerkerneval(cgrph, kern, sigma, [x(:) y(:)].', opts);
                u = reshape(u, size(x));
                u = real(u);
            end
            u.bulk = @u_bulk;

        end

    end

end

function A = edgewise(cgrph, f)
    blocks = arrayfun(f, cgrph.echnks, 'UniformOutput', false);
    A = matlab.internal.math.blkdiag(blocks{:});
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

function afuns = arclengthfun_graph(cgrph)
    vec = @(x) x(:);
    afuns = arrayfun(@(x) vec(arclengthfun(x)), cgrph.echnks, 'UniformOutput', false);
    afuns = afuns.';
end
