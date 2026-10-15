classdef TestPerfusionOp < matlab.unittest.TestCase
%TESTPERFUSIONOP   Unit tests for PERFUSIONOP on small synthetic networks.
%   The networks are a single straight segment and a Y-shaped tree with
%   its root at vertex 1. Tests check the boundary and vertex conditions of
%   the solution, flux conservation, the perfectly conducting limit
%   (xi = 0), agreement between solvers, and the arc-length integration
%   helpers, with xi and λ given as scalars, per-edge constants, smooth
%   functions of arc length, and nodal values.

    properties
        segment  % Single straight edge from vertex 1 (root) to vertex 2
        tree     % Y-shaped tree: root -> junction -> two tips
    end

    properties (Constant)
        params = {'xi', 1, 'phi', 2, 'P', 0.03}
        tol = 1e-10
        % Smooth, positive xi(s) and λ(s) on each edge of the tree
        xif = {@(s) 1 + 0.5*sin(3*s), @(s) 0.5 + s.^2, @(s) 2 - 0.8*cos(2*s)}
        lamf = {@(s) 0.1 + 0.05*s, @(s) 0.2*ones(size(s)), @(s) 0.05 + 0.1*s.^2}
    end

    methods (TestClassSetup)
        function buildNetworks(tc)
            cparams = struct('nover', 2);
            % Edges are given as [start; end] vertex indices
            verts = [0 0; 1 0].';
            tc.segment = chunkgraph(verts, [1; 2], [], cparams);

            verts = [0 0; 0 1; -0.7 1.7; 0.7 1.7].';
            tc.tree = chunkgraph(verts, [1 2 2; 2 3 4], [], cparams);
        end
    end

    methods (Test)

        function sizeMatchesUnknowns(tc)
            op = perfusionop(tc.tree, tc.params{:});
            tc.verifyEqual(size(op), [op.ndof op.ndof]);
            tc.verifyEqual(op.ndof, op.npt + 2*op.nedges);
            opo = perfusionop(tc.tree, tc.params{:}, outer=true);
            tc.verifyEqual(opo.ndof, opo.npt + 2*opo.nedges + opo.outer.npt);
        end

        function rootValueIsOne(tc)
            for alpha = [0 0.1]
                op = perfusionop(tc.tree, tc.params{:}, alpha=alpha);
                u = op.solve(method='direct');
                v = edgeEndValues(op, u, 1);
                tc.verifyEqual(v(1), 1, AbsTol=tc.tol);
            end
        end

        function networkIsContinuousAtJunction(tc)
            op = perfusionop(tc.tree, tc.params{:}, alpha=0.1);
            u = op.solve(method='direct');
            v1 = edgeEndValues(op, u, 1);  % Ends at the junction
            v2 = edgeEndValues(op, u, 2);  % Starts at the junction
            v3 = edgeEndValues(op, u, 3);  % Starts at the junction
            tc.verifyEqual(v2(1), v1(2), AbsTol=tc.tol);
            tc.verifyEqual(v3(1), v1(2), AbsTol=tc.tol);
        end

        function concentrationDecreasesFromRoot(tc)
            op = perfusionop(tc.segment, tc.params{:}, alpha=0.1);
            u = op.solve(method='direct');
            v = u.channel{1};
            tc.verifyTrue(all(v > 0 & v < 1));
            tc.verifyTrue(all(diff(v) < 0));
        end

        function ingressFluxEqualsTissueFlux(tc)
            % Integrating the conservation law along the network and using
            % the Kirchhoff and tip conditions, the flux entering at the
            % root, -a v'(0) = -b_1, equals the total flux into the
            % tissue, the integral of σ.
            for alpha = [0 0.1]
                op = perfusionop(tc.tree, tc.params{:}, alpha=alpha);
                [~, vars] = op.solve(method='direct');
                ingress = -vars.b(1);
                tissue = sum(vars.sigma .* op.cgrph.wts(:));
                tc.verifyEqual(tissue, ingress, RelTol=1e-10);
                tc.verifyGreaterThan(ingress, 0);
            end
        end

        function conservationWithAnyRoots(tc)
            % Roots at the start of an edge (1), at the end of an edge (4),
            % at the degree-3 junction (2), and at two vertices (1 and 3).
            % v = 1 on every edge at every root, and the ingress flux summed
            % over all roots and incident edges equals ∫σ and the absorption.
            for root = {1, 4, 2, [1 3]}
                r = root{1};
                op = perfusionop(tc.tree, tc.params{:}, alpha=0.1, root=r);
                [u, vars] = op.solve(method='direct');
                for vert = r
                    edge_ids = op.cgrph.vstruc{vert}{1};
                    edge_dir = op.cgrph.vstruc{vert}{2};
                    for k = 1:numel(edge_ids)
                        v = edgeEndValues(op, u, edge_ids(k));
                        tc.verifyEqual(v((edge_dir(k)+3)/2), 1, ...
                            sprintf('root=%s edge %d', mat2str(r), edge_ids(k)), AbsTol=tc.tol);
                    end
                end
                vals = [];
                evalc('vals = op.check_conservation(u, vars, ngrid=200);');
                msg = sprintf('root=%s', mat2str(r));
                tc.verifyEqual(vals.ingress, vals.int_sigma, msg, RelTol=1e-10);
                tc.verifyEqual(vals.ingress_dv, vals.int_sigma, msg, RelTol=1e-6);
                tc.verifyEqual(vals.absorption, vals.int_sigma, msg, RelTol=1e-2);
            end
        end

        function absorptionEqualsIngressWithOuterBoundary(tc)
            % With a no-flux outer boundary, all O2 entering at the root is
            % absorbed in the disk: φ² ∫ u = ingress flux. The area
            % integral uses a polar grid, so the tolerance is loose.
            op = perfusionop(tc.tree, tc.params{:}, alpha=0.1, outer=true);
            [u, vars] = op.solve(method='direct');
            ingress = -vars.b(1);

            r = op.outer.r(:,:);
            c = mean(r, 2);
            R = mean(vecnorm(r - c));
            [rr, wr] = lege.exps(200);
            rr = R*(rr+1)/2;  wr = R*wr/2;
            nth = 400;
            th = 2*pi*(0:nth-1)/nth;
            [RR, TH] = ndgrid(rr, th);
            uu = u.bulk(c(1) + RR.*cos(TH), c(2) + RR.*sin(TH), ...
                struct(accel=true, forcesmooth=true));
            absorbed = op.phi^2 * (2*pi/nth) * sum(wr(:) .* rr(:) .* sum(uu, 2));
            tc.verifyEqual(absorbed, ingress, RelTol=1e-3);
        end

        function perfectlyConductingNetwork(tc)
            % With xi = 0 (infinite channel area), v = 1 on the whole
            % network and σ solves the first-kind equation λσ + Sσ = η.
            % Compare with a direct solve of that equation.
            phi = 2;  P = 0.03;
            S = chunkermat(tc.tree, kernel('helm', 's', 1i*phi), struct(rcip=false));
            for alpha = [0 0.1]
                op = perfusionop(tc.tree, xi=0, phi=phi, P=P, alpha=alpha);
                [u, vars] = op.solve(method='direct');
                tc.verifyEqual(vertcat(u.channel{:}), ones(op.npt, 1), AbsTol=tc.tol);
                sigma = (S + alpha*eye(op.npt)) \ (P*ones(op.npt, 1));
                tc.verifyLessThan(norm(vars.sigma - sigma) / norm(sigma), 1e-12);
                [~, vars2] = op.solve(method='gmres', schur=true, precond=(alpha == 0), tol=1e-12);
                tc.verifyLessThan(norm(vars2.sigma - vars.sigma) / norm(vars.sigma), 1e-9);
            end
        end

        function solversAgree(tc)
            for alpha = [0 0.1]
                for outer = [false true]
                    op = perfusionop(tc.tree, tc.params{:}, alpha=alpha, outer=outer);
                    precond = (alpha == 0);
                    [~, ref] = op.solve(method='direct');
                    configs = { ...
                        {'method', 'direct', 'schur', true}, ...
                        {'method', 'gmres',  'schur', false, 'precond', precond}, ...
                        {'method', 'gmres',  'schur', true,  'precond', precond}};
                    for k = 1:numel(configs)
                        [~, vars] = op.solve(configs{k}{:}, tol=1e-12);
                        err = norm(vars.sigma - ref.sigma) / norm(ref.sigma);
                        tc.verifyLessThan(err, 1e-9, ...
                            sprintf('alpha=%g outer=%d %s', alpha, outer, strjoin(string(configs{k}), ' ')));
                    end
                end
            end
        end

        function fmmSolversAgree(tc)
            % Exercises the FMM code paths, including the outer-boundary
            % preconditioner without a Schur complement.
            tc.assumeTrue(exist('hfmm2d', 'file') > 0, 'fmm2d is not available.');
            for outer = [false true]
                op = perfusionop(tc.tree, tc.params{:}, alpha=0.1, outer=outer);
                [~, ref] = op.solve(method='direct');
                for schur = [false true]
                    [~, vars] = op.solve(method='gmres', fmm=true, schur=schur, tol=1e-12);
                    err = norm(vars.sigma - ref.sigma) / norm(ref.sigma);
                    tc.verifyLessThan(err, 1e-8, sprintf('outer=%d schur=%d', outer, schur));
                end
            end
        end

        function gmresReturnsInfo(tc)
            op = perfusionop(tc.segment, tc.params{:}, alpha=0.1);
            [~, ~, info] = op.solve(method='gmres', tol=1e-10);
            tc.verifyEqual(info.flag, 0);
            tc.verifyLessThan(info.relres, 1e-10);
            tc.verifyNotEmpty(info.resvec);
        end

        function antiderivativeOfConstant(tc)
            % Q 1 = s and Q² 1 = s²/2 on every edge.
            op = perfusionop(tc.tree, tc.params{:});
            one = ones(op.npt, 1);
            s = vertcat(op.s{:});
            tc.verifyEqual(op.iintegral(op.cgrph, one, 1), s,      AbsTol=1e-13);
            tc.verifyEqual(op.iintegral(op.cgrph, one, 2), s.^2/2, AbsTol=1e-13);
        end

        function integrationRows(tc)
            % With xi(s) = 1 + s on each edge: introw1{k}*f = ∫ f,
            % Qxi{k} = s + s²/2, and introwxi{k}*f = (Q xi Q f)(L).
            xi = repmat({@(s) 1 + s}, 1, 3);
            op = perfusionop(tc.tree, xi=xi, phi=2, P=0.03);
            for k = 1:op.nedges
                s = op.s{k};
                f = cos(3*s);
                L = op.L(k);
                QxiQf = (1 - cos(3*L))/9 + (-L*cos(3*L)/3 + sin(3*L)/9)/3;
                tc.verifyEqual(op.introw1{k}*f, sin(3*L)/3, AbsTol=1e-12);
                tc.verifyEqual(op.Qxi{k}, s + s.^2/2, AbsTol=1e-12);
                tc.verifyEqual(op.introwxi{k}*f, QxiQf, AbsTol=1e-12);
            end
        end

        function coefficientFormsAgree(tc)
            % A scalar, one value per edge, one function per edge, and one
            % value per node all give the same operator.
            ref = perfusionop(tc.tree, xi=0.7, alpha=0.1, phi=2, P=0.03);
            [~, vref] = ref.solve(method='direct');
            npt = ref.npt;
            xis = {0.7*[1 1 1], {@(s) 0.7 + 0*s, 0.7, @(s) 0.7*ones(size(s))}, 0.7*ones(npt, 1)};
            alphas = {0.1*[1 1 1], {0.1, @(s) 0.1 + 0*s, 0.1}, 0.1*ones(npt, 1)};
            for k = 1:numel(xis)
                op = perfusionop(tc.tree, xi=xis{k}, alpha=alphas{k}, phi=2, P=0.03);
                tc.verifyEqual(op.xi, ref.xi);
                tc.verifyEqual(op.alpha, ref.alpha);
                [~, vars] = op.solve(method='direct');
                tc.verifyLessThan(norm(vars.sigma - vref.sigma) / norm(vref.sigma), 1e-13);
            end
        end

        function perEdgeCoefficients(tc)
            % Different constants on each edge: v = 1 at the root, v is
            % continuous at the junction, the flux Qσ + b balances there and
            % vanishes at the tips, and the ingress flux equals ∫σ.
            op = perfusionop(tc.tree, xi=[1 0.3 3], alpha=[0.1 0.05 0.2], phi=2, P=0.03);
            [u, vars] = op.solve(method='direct');
            v1 = edgeEndValues(op, u, 1);
            v2 = edgeEndValues(op, u, 2);
            v3 = edgeEndValues(op, u, 3);
            tc.verifyEqual(v1(1), 1, AbsTol=tc.tol);
            tc.verifyEqual(v2(1), v1(2), AbsTol=tc.tol);
            tc.verifyEqual(v3(1), v1(2), AbsTol=tc.tol);
            intsigma = arrayfun(@(k) op.introw1{k} * ...
                vars.sigma(op.sigma_starts(k):op.sigma_starts(k+1)-1), 1:3);
            fluxend = intsigma + vars.b.';
            tc.verifyEqual(fluxend(1) - vars.b(2) - vars.b(3), 0, AbsTol=tc.tol);
            tc.verifyEqual(fluxend(2:3), [0 0], AbsTol=tc.tol);
            tc.verifyEqual(-vars.b(1), sum(intsigma), RelTol=1e-10);
        end

        function networkEquationWithVariableCoefficients(tc)
            % With smooth xi(s) and λ(s), a v' = v'/xi equals Qσ + b at the
            % nodes, where v' is computed by differentiating v on each chunk.
            op = perfusionop(tc.tree, xi=tc.xif, alpha=tc.lamf, phi=2, P=0.03);
            [u, vars] = op.solve(method='direct');
            n = op.n;
            D = lege.dermat(n);
            for k = 1:op.nedges
                chnkr = op.cgrph.echnks(k);
                idx = op.sigma_starts(k):op.sigma_starts(k+1)-1;
                ds = reshape(vecnorm(chnkr.d(:,:)), n, []);
                dv = (D * reshape(u.channel{k}, n, [])) ./ ds;
                flux = op.iintegral(chnkr, vars.sigma(idx), 1) + vars.b(k);
                tc.verifyLessThan(norm(dv(:)./op.xi(idx) - flux, inf) / norm(flux, inf), 1e-6, ...
                    sprintf('edge %d', k));
            end
        end

        function variableCoefficientSolversAgree(tc)
            for outer = [false true]
                op = perfusionop(tc.tree, xi=tc.xif, alpha=tc.lamf, phi=2, P=0.03, outer=outer);
                [~, ref] = op.solve(method='direct');
                configs = {{'method', 'direct', 'schur', true}, ...
                           {'method', 'gmres', 'schur', false}, ...
                           {'method', 'gmres', 'schur', true}};
                if ( exist('hfmm2d', 'file') > 0 )
                    configs{end+1} = {'method', 'gmres', 'schur', true, 'fmm', true};
                end
                for k = 1:numel(configs)
                    [~, vars] = op.solve(configs{k}{:}, tol=1e-12);
                    err = norm(vars.sigma - ref.sigma) / norm(ref.sigma);
                    tc.verifyLessThan(err, 1e-8, ...
                        sprintf('outer=%d %s', outer, strjoin(string(configs{k}), ' ')));
                end
            end
        end

        function convergenceWithVariableCoefficients(tc)
            % Refining dyadically into the vertices, u at a target converges
            % to high accuracy with smooth xi(s) and λ(s).
            cgrph = tc.tree;
            targ = [0.35 0.8];
            nlev = 7;
            vals = zeros(nlev, 1);
            for lev = 1:nlev
                op = perfusionop(cgrph, xi=tc.xif, alpha=tc.lamf, phi=2, P=0.03);
                u = op.solve(method='direct');
                vals(lev) = u.bulk(targ(1), targ(2));
                cgrph = refineIntoVertices(cgrph);
            end
            diffs = abs(diff(vals));
            tc.verifyTrue(all(diff(diffs) < 0), 'differences should decrease');
            tc.verifyLessThan(diffs(end), 1e-11);
        end

        function perfectlyConductingEdge(tc)
            % xi = 0 on edge 2 only: v is constant there and equal to v at
            % the junction, and the ingress flux still equals ∫σ.
            op = perfusionop(tc.tree, xi=[1 0 1], alpha=0.1, phi=2, P=0.03);
            [u, vars] = op.solve(method='direct');
            v1 = edgeEndValues(op, u, 1);
            tc.verifyEqual(u.channel{2}, v1(2)*ones(size(u.channel{2})), AbsTol=tc.tol);
            tc.verifyEqual(-vars.b(1), sum(vars.sigma .* op.cgrph.wts(:)), RelTol=1e-10);
        end

        function badCoefficients(tc)
            id = 'perfusionop:badCoefficient';
            make = @(varargin) perfusionop(tc.tree, varargin{:});
            tc.verifyError(@() make(xi=[1 2]), id);                 % wrong length
            tc.verifyError(@() make(xi=-1), id);                    % negative
            tc.verifyError(@() make(alpha=NaN), id);                % not finite
            tc.verifyError(@() make(alpha=[0.1 -0.1 0.1]), id);     % negative on an edge
            tc.verifyError(@() make(xi={1, 2}), id);                % too few cells
            tc.verifyError(@() make(xi={@(s) [1 2], 1, 1}), id);    % wrong number of values
            tc.verifyError(@() make(xi={'a', 1, 1}), id);           % not a function or scalar
            tc.verifyError(@() make(xi='a'), id);                   % not numeric
        end

    end

end

function v = edgeEndValues(op, u, k)
%EDGEENDVALUES   Network concentration at the start and end of edge K.
%   Interpolates U.channel{K} from the first and last chunks of the edge to
%   their endpoints.
    chnkr = op.cgrph.echnks(k);
    n = chnkr.k;
    vals = reshape(u.channel{k}, n, []);
    first = find(chnkr.adj(1,:) <= 0, 1);
    last  = find(chnkr.adj(2,:) <= 0, 1);
    v = [lege.matrin(n, -1) * vals(:,first), ...
         lege.matrin(n,  1) * vals(:,last)];
end

function cgrph = refineIntoVertices(cgrph)
%REFINEINTOVERTICES   Split the first and last chunk of every edge.
%   Halves the chunks touching each vertex, where σ is singular; the
%   chunkgraph refinement then balances neighboring chunks.
    split = cell(1, numel(cgrph.echnks));
    for k = 1:numel(cgrph.echnks)
        adj = cgrph.echnks(k).adj;
        split{k} = unique([find(adj(1,:) <= 0, 1), find(adj(2,:) <= 0, 1)]);
    end
    cgrph = refine(cgrph, struct('splitchunks', {split}));
end
