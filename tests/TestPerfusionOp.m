classdef TestPerfusionOp < matlab.unittest.TestCase
%TESTPERFUSIONOP   Unit tests for PERFUSIONOP on small synthetic networks.
%   The networks are a single straight segment and a Y-shaped tree with
%   its root at vertex 1. Tests check the boundary and vertex conditions of
%   the solution, flux conservation, agreement between solvers, and the
%   arc-length integration helpers.

    properties
        segment  % Single straight edge from vertex 1 (root) to vertex 2
        tree     % Y-shaped tree: root -> junction -> two tips
    end

    properties (Constant)
        params = {'xi', 1, 'phi', 2, 'P', 0.03}
        tol = 1e-10
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
            % root, -a v'(0) = -b_1/xi, equals the total flux into the
            % tissue, the integral of σ.
            for alpha = [0 0.1]
                op = perfusionop(tc.tree, tc.params{:}, alpha=alpha);
                [~, vars] = op.solve(method='direct');
                ingress = -vars.b(1) / op.xi;
                tissue = sum(vars.sigma .* op.cgrph.wts(:));
                tc.verifyEqual(tissue, ingress, RelTol=1e-10);
                tc.verifyGreaterThan(ingress, 0);
            end
        end

        function absorptionEqualsIngressWithOuterBoundary(tc)
            % With a no-flux outer boundary, all O2 entering at the root is
            % absorbed in the disk: φ² ∫ u = ingress flux. The area
            % integral uses a polar grid, so the tolerance is loose.
            op = perfusionop(tc.tree, tc.params{:}, alpha=0.1, outer=true);
            [u, vars] = op.solve(method='direct');
            ingress = -vars.b(1) / op.xi;

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
            % introw1{k}*f = ∫ f and introw2{k}*f = (Q² f)(L) on edge k.
            op = perfusionop(tc.tree, tc.params{:});
            for k = 1:op.nedges
                s = op.s{k};
                f = cos(3*s);
                L = op.L(k);
                tc.verifyEqual(op.introw1{k}*f, sin(3*L)/3, AbsTol=1e-12);
                tc.verifyEqual(op.introw2{k}*f, (1 - cos(3*L))/9, AbsTol=1e-12);
            end
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
