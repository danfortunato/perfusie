classdef perfusionop %#ok<*PROP,*PROPLC>
%PERFUSIONOP   Integral-equation solver for network perfusion of thin tissue.
%   OP = PERFUSIONOP(CGRPH) sets up the coupled PDE-ODE perfusion model on
%   the network CGRPH, a chunkgraph embedded in the plane:
%
%     (Δ - φ²) u = 0                in the tissue, off the network,
%     λ [u_n] - (η v - u) = 0       on each edge (membrane permeation),
%     (a v')' - [u_n] = 0           on each edge (conservation),
%     v continuous, Σ ±a v' = 0     at each non-root vertex (Kirchhoff),
%     v = 1                         at each root vertex,
%
%   where u is the tissue concentration, v the network concentration,
%   [u_n] the jump in normal derivative across an edge, and a = 1/xi the
%   scaled channel cross-sectional area. With OUTER = true, the tissue is
%   bounded by a circle with zero-flux conditions; otherwise it fills the
%   plane.
%
%   The tissue concentration is represented as a single-layer potential
%   u = S σ on the network (plus S τ on the outer boundary), and on each
%   edge v = xi Q² σ + a_k + b_k s, where Q is arc-length integration from
%   the start of the edge. This yields a second-kind-like system for σ (and
%   τ) coupled to two constants (a_k, b_k) per edge.
%
%   OP = PERFUSIONOP(CGRPH, Name=Value) sets parameters:
%     xi    - inverse scaled channel area, 1/a (default 0).
%     alpha - inverse membrane permeability, λ (default 0).
%     phi   - Yukawa parameter φ; 1/φ is the decay length in the tissue
%             (default 0).
%     P     - tissue-air partition (Henry's) constant, η (default 0.03).
%     root  - indices of the root (ingress) vertices of CGRPH (default 1).
%     outer - if true, enclose the network in a no-flux circle 1.3 times
%             the radius of its minimal bounding circle (default false).
%
%   Example:
%     cgrph = graphtools.fit(graphtools.gallery('L1_1_Tr9L'));
%     op = perfusionop(cgrph, xi=0.004, phi=0.1);
%     [u, vars] = op.solve();
%     plot_perfusion(op, u)
%
%   PERFUSIONOP methods:
%     solve              - Solve the model.
%     discretize         - Assemble the discretized system.
%     check_conservation - Print a check of global O2 conservation.
%     size               - Number of unknowns.
%
%   See also PLOT_PERFUSION, GRAPHTOOLS.FIT.

    properties
        cgrph        % Network, a chunkgraph
        root         % Indices of the root vertices
        xi           % Inverse scaled channel area
        alpha        % Inverse membrane permeability
        phi          % Yukawa parameter
        P            % Henry's constant
        nverts       % Number of vertices
        nedges       % Number of edges
        n            % Legendre order of each chunk
        nch          % Total number of chunks on the network
        npt          % Total number of discretization nodes on the network
        nbc          % Number of edge constants (two per edge)
        ndof         % Total number of unknowns
        sigma_starts % Index of the first node of each edge in σ (plus end)
        L            % Arc length of each edge
        s            % Arc-length coordinate at the nodes of each edge
        introw1      % Row vectors giving ∫σ over each edge
        introw2      % Row vectors giving (Q²σ)(L) on each edge
        iintegral    % Handle to the arc-length antiderivative routine
        outer        % Outer boundary chunker, or [] if none
    end

    methods

        function op = perfusionop(cgrph, opts)
        %PERFUSIONOP   Construct a perfusion operator. See PERFUSIONOP.
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
        %SIZE   Size of the discretized system.
        %   SZ = SIZE(OP) returns [NDOF NDOF], where NDOF is the number of
        %   unknowns. [M, N] = SIZE(OP) and SIZE(OP, DIM) work as for
        %   matrices.
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
        %DISCRETIZE   Assemble the discretized system.
        %   [A, RHS, R] = DISCRETIZE(OP) returns the dense system matrix A,
        %   right-hand side RHS, and right preconditioner R, with unknowns
        %   ordered [σ; a; b; τ] (τ only if OP has an outer boundary). The
        %   system to solve is A*R*y = RHS with x = R*y.
        %
        %   DISCRETIZE(OP, Name=Value) sets options:
        %     blocks  - if true, return A, R, and RHS as cell arrays of
        %               blocks, with block rows/columns for (σ, [a; b], τ)
        %               (default false).
        %     fmm     - if true, return the dense-operator blocks as
        %               function handles f(x, tol) applied with the FMM.
        %               Implies blocks = true (default false).
        %     precond - if true, right-precondition the σ block with the
        %               hypersingular operator D' on each edge; useful when
        %               alpha is small (default false).
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

                    % Target the merged edge chunkers (same points, same order)
                    % since chunkerkernevalmat mishandles chunkgraph targets.
                    cormat = chunkerkernevalmat(op.outer, Skern, echnks, struct(corrections=true));
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
                    % Target the merged edge chunkers (same points, same order)
                    % since chunkerkernevalmat mishandles chunkgraph targets.
                    out2in = chunkerkernevalmat(op.outer, Skern,  merge(op.cgrph.echnks));
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

        function [u, vars, info] = solve(op, opts)
        %SOLVE   Solve the perfusion model.
        %   [U, VARS] = SOLVE(OP) solves the model and returns the solution
        %   in a struct U with fields
        %     channel - cell array with the network concentration v at the
        %               nodes of each edge,
        %     bulk    - function handle U.bulk(X, Y, OPTS) that evaluates
        %               the tissue concentration u at points (X, Y), passing
        %               OPTS to chunkerkerneval,
        %   and the raw unknowns in a struct VARS with fields sigma, a, b,
        %   and tau (if OP has an outer boundary).
        %
        %   [U, VARS, INFO] = SOLVE(OP, ...) also returns the GMRES flag,
        %   relres, iter, and resvec in the struct INFO (empty for direct
        %   solves).
        %
        %   SOLVE(OP, Name=Value) sets options:
        %     method   - 'direct', 'gmres', or 'auto' (default). 'auto'
        %                uses a direct solve below 10000 unknowns and GMRES
        %                with the FMM and Schur complement above.
        %     fmm      - if true, apply dense blocks with the FMM; GMRES
        %                only (default false).
        %     schur    - if true, eliminate the edge constants (a, b) and
        %                solve the Schur complement for σ (and τ)
        %                (default false).
        %     precond  - if true, use D' right preconditioning; see
        %                DISCRETIZE (default false).
        %     showcond - if true, print the condition number of the
        %                (preconditioned) system. Disables the FMM
        %                (default false).
        %     tol      - GMRES tolerance (default 1e-10).
        %     maxit    - maximum GMRES iterations (default OP.ndof).

            arguments
                op
                opts.method   = 'auto'
                opts.fmm      = false
                opts.schur    = false
                opts.precond  = false
                opts.showcond = false
                opts.tol      = 1e-10
                opts.maxit    = op.ndof
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

            if ( opts.showcond )
                opts.fmm = false;
            end

            info = struct();
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
                        R = @(x) [ R{1,1}(x(ii),tol) ; R{2,2}*x(bb) ; R{3,3}(x(jj),tol) ];
                        rhs = [rhs{1} ; rhs{2} ; rhs{3}];
                    else
                        A = @(x) [ A{1,1}(x(ii),tol) + A{1,2}*x(bb) ;
                                   A{2,1}*x(ii)      + A{2,2}*x(bb) ];
                        R = @(x) [ R{1,1}(x(ii),tol) ; R{2,2}*x(bb) ];
                        rhs = [rhs{1} ; rhs{2}];
                    end
                end
                if ( opts.showcond )
                    fprintf('cond(A) = %g (%g)\n', cond(full(A*R)), skeel(full(A*R), 2));
                end
                if strcmpi(opts.method, 'gmres')
                    [sigma_a_b, flag, relres, iter, resvec] = solve(A, rhs, R);
                    info.flag = flag;
                    info.relres = relres;
                    info.iter = iter;
                    info.resvec = resvec;
                else
                    sigma_a_b = solve(A, rhs, R);
                end
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
                if ( opts.showcond )
                    fprintf('cond(S) = %g (%g)\n', cond(full(Sii*Rii)), skeel(full(Sii*Rii), 2));
                    % fprintf('cond(Abb) = %g (%g)\n', cond(full(A{2,2})), skeel(full(A{2,2})));
                end
                if strcmpi(opts.method, 'gmres')
                    [sigma, flag, relres, iter, resvec] = solve(Sii, fi, Rii);
                    info.flag = flag;
                    info.relres = relres;
                    info.iter = iter;
                    info.resvec = resvec;
                else
                    sigma = solve(Sii, fi, Rii);
                end
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
            %U_BULK   Evaluate the tissue concentration at points (X, Y).
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
        %CHECK_CONSERVATION   Print a check of global O2 conservation.
        %   CHECK_CONSERVATION(OP, U, VARS) prints four quantities that
        %   should agree: the ingress flux at the root computed (1) from the
        %   edge constant b and (2) by differentiating v, (3) the total flux
        %   ∫σ into the tissue, and (4) the total absorption φ²∫u, computed
        %   with a uniform grid over the tissue. It then prints the relative
        %   difference between (1) and (4).
        %
        %   CHECK_CONSERVATION(OP, U, VARS, ngrid=N) uses an N-by-N grid for
        %   (4) (default 800).
        %
        %   This assumes a single root vertex, incident to edge 1 at the
        %   start of that edge.

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
%MYINTMAT   Dense matrix of arc-length integration along a chunker.
%   IMAT = MYINTMAT(CHNKR) returns the matrix that maps values of f at the
%   nodes of CHNKR to values of its arc-length antiderivative Q f at the
%   same nodes. The constant of integration is chosen so that Q f is zero
%   at the left end of the first chunk. Chunks are visited in the order
%   given by CHNKR.adj.

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
%DINTEGRAL   Definite arc-length integral over a chunker.
%   INT = DINTEGRAL(CHNKR, F) returns the integral of F, given at the nodes
%   of CHNKR, over the whole chunker.
    int = sum(chnkr.wts(:) .* f);
end

function int = iintegral(chnkr, f, p)
%IINTEGRAL   Repeated arc-length antiderivative along each edge.
%   INT = IINTEGRAL(CHNKR, F) returns Q F at the nodes of CHNKR, where
%   (Q F)(t) is the integral of F from the start of the chunker to t.
%   INT = IINTEGRAL(CHNKR, F, P) applies Q P times (default 1).
%
%   If CHNKR is a chunkgraph, Q is applied separately on each edge, with F
%   ordered edge by edge as in the chunkgraph. Chunks within each chunker
%   are assumed to be stored in order along the curve.

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
%IINTEGRALT   Transpose of the discrete arc-length antiderivative.
%   INT = IINTEGRALT(CHNKR, F, P) applies the transpose of the discretized
%   operator Q from IINTEGRAL, P times (default 1). Q^T maps quadrature-
%   weighted values to the integral from each node to the end of the curve.
%   Accepts a chunker or a chunkgraph, as for IINTEGRAL.

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
%INTROW   Row vector evaluating a repeated antiderivative at the end.
%   R = INTROW(CHNKR, P) returns the row vector R such that R*F equals
%   (Q^P F)(L), the P-fold arc-length antiderivative of F evaluated at the
%   end of CHNKR (default P = 1, so that R*F is the integral of F).

    if ( nargin < 2 )
        p = 1;
    end

    w = chnkr.wts(:);
    introw = iintegralT(chnkr, w, p-1).';

end

function [x, flag, relres, iter, resvec] = right_gmres(A, b, restart, tol, maxit, P)
%RIGHT_GMRES   GMRES with right preconditioning.
%   [X, FLAG, RELRES, ITER, RESVEC] = RIGHT_GMRES(A, B, RESTART, TOL,
%   MAXIT, P) solves A*X = B by running GMRES on A*P*Y = B and returning
%   X = P*Y. A and P may be matrices or function handles. The other inputs
%   and outputs are as for GMRES.

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
%BLOCKDIAG   Sparse block-diagonal part of a matrix.
%   B = BLOCKDIAG(A, N) returns a sparse matrix holding the N-by-N diagonal
%   blocks of the square matrix A, whose size must be a multiple of N.

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
