classdef TestGraphTools < matlab.unittest.TestCase
%TESTGRAPHTOOLS   Unit tests for the GRAPHTOOLS package on a bundled cell.

    properties
        G           % Imported graph of a small L1 cell
        Gs          % Smoothed graph
        superedges  % Super-edges of Gs
    end

    properties (Constant)
        name = 'L1_1_Tr9L'
    end

    methods (TestClassSetup)
        function loadCell(tc)
            tc.G = graphtools.gallery(tc.name);
            [tc.Gs, tc.superedges] = graphtools.smooth_graph(tc.G, 1);
        end
    end

    methods (Test)

        function galleryReturnsTree(tc)
            G = tc.G;
            n = numnodes(G);
            tc.verifyGreaterThan(n, 1);
            tc.verifySize(G.Nodes.coords, [n 3]);
            tc.verifySize(G.Nodes.radius, [n 1]);
            tc.verifyEqual(numedges(G), n - 1);
            tc.verifyEqual(max(conncomp(G)), 1);
            % A few edges have zero length, where a traced path starts on an
            % existing point.
            tc.verifyGreaterThanOrEqual(G.Edges.length, 0);
        end

        function galleryByIndex(tc)
            files = dir(fullfile(dataDir(), '*.traces'));
            G1 = graphtools.gallery(1);
            Gref = graphtools.import_graph(fullfile(dataDir(), files(1).name));
            tc.verifyEqual(G1.Nodes.coords, Gref.Nodes.coords);
        end

        function galleryUnknownName(tc)
            tc.verifyError(@() graphtools.gallery('no_such_cell'), ...
                'graphtools:gallery:notFound');
        end

        function everyTraceHasAnImage(tc)
            % Every trace has a matching PDF, except L1_16 (see data/README.md).
            traces = erase(string({dir(fullfile(dataDir(), '*.traces')).name}), '.traces');
            pdfs = erase(string({dir(fullfile(dataDir(), '*.pdf')).name}), '.pdf');
            tc.verifyEqual(sort(setdiff(traces, pdfs)), ["L1_16_Tr9L" "L1_16_Tr9R"]);
            tc.verifyEmpty(setdiff(pdfs, traces));
        end

        function smoothingKeepsBranchPointsAndTips(tc)
            anchors = degree(tc.G) ~= 2;
            tc.verifyEqual(tc.Gs.Nodes.coords(anchors,1:2), tc.G.Nodes.coords(anchors,1:2));
            tc.verifyEqual(numnodes(tc.Gs), numnodes(tc.G));
            tc.verifyEqual(tc.Gs.Edges.EndNodes, tc.G.Edges.EndNodes);
        end

        function superedgesPartitionEdges(tc)
            % In a tree, the super-edges cover each edge exactly once and
            % run between nodes of degree ~= 2.
            se = tc.superedges;
            deg = degree(tc.G);
            tc.verifyEqual(sum(arrayfun(@(e) numel(e.nodes) - 1, se)), numedges(tc.G));
            tc.verifyFalse(any([se.is_cycle]));
            ends = cell2mat(arrayfun(@(e) e.nodes([1 end]), se, UniformOutput=false));
            tc.verifyTrue(all(deg(ends) ~= 2, 'all'));
        end

        function superedgeGraphMatchesSuperedges(tc)
            S = graphtools.superedge_graph(tc.Gs, tc.superedges, directed=true);
            tc.verifyClass(S, 'digraph');
            tc.verifyEqual(numedges(S), numel(tc.superedges));
            tc.verifyEqual(sort(S.Edges.superedge_id), (1:numel(tc.superedges)).');
            tc.verifyEqual(S.Nodes.orig_id(1), 1);  % The root comes first
        end

        function fitEndpointsMatchVertices(tc)
            cgrph = graphtools.graph_to_chunkgraph(tc.Gs, tc.superedges, tol=1e-3);
            tc.verifyEqual(numel(cgrph.echnks), numel(tc.superedges));
            for k = 1:numel(cgrph.echnks)
                chnkr = cgrph.echnks(k);
                ends = chunkends(chnkr, [1 chnkr.nch]);
                xy = tc.Gs.Nodes.coords(tc.superedges(k).nodes([1 end]), 1:2).';
                tc.verifyEqual(ends(:,[1 4]), xy, sprintf('edge %d', k), AbsTol=1e-8);
            end
        end

        function fitPutsRootAtVertexOne(tc)
            % PERFUSIONOP takes root=1 by default, which relies on this.
            cgrph = graphtools.fit(tc.G, tol=1e-3);
            tc.verifyEqual(cgrph.verts(:,1), tc.G.Nodes.coords(1,1:2).', AbsTol=1e-12);
        end

        function fitAcceptsFileOrGraph(tc)
            file = fullfile(dataDir(), [tc.name '.traces']);
            c1 = graphtools.fit(file, tol=1e-3);
            c2 = graphtools.fit(tc.G, tol=1e-3);
            tc.verifyEqual(c1.r, c2.r);
        end

        function fitOversampling(tc)
            c0 = graphtools.fit(tc.G, tol=1e-3, nover=0);
            c1 = graphtools.fit(tc.G, tol=1e-3, nover=1);
            tc.verifyGreaterThan(c1.npt, c0.npt);
            tc.verifyEqual(c1.verts, c0.verts);
        end

        function coarsenRemovesDegreeTwoNodes(tc)
            S = graphtools.coarsen_graph(tc.G);
            deg = degree(S);
            root = find(S.Nodes.orig_id == 1);
            tc.verifyEqual(root, 1);
            tc.verifyFalse(any(deg(2:end) == 2));
            tc.verifyEqual(nnz(degree(tc.G) ~= 2), numnodes(S) - (degree(tc.G, 1) == 2));
            tc.verifyEqual(S.Nodes.coords, tc.G.Nodes.coords(S.Nodes.orig_id,:));
        end

        function rasterizeCoversNodes(tc)
            dx = 0.5;
            [BW, xg, yg] = graphtools.rasterize_graph(tc.G, dx);
            tc.verifySize(BW, [numel(yg) numel(xg)]);
            xy = tc.G.Nodes.coords(:,1:2);
            cols = round((xy(:,1) - xg(1))/dx) + 1;
            rows = round((xy(:,2) - yg(1))/dx) + 1;
            tc.verifyTrue(all(BW(sub2ind(size(BW), rows, cols))));
        end

        function digraphToGraph(tc)
            D = digraph([1 2 3], [2 3 1]);
            G = graphtools.digraph_to_graph(D);
            tc.verifyClass(G, 'graph');
            tc.verifyEqual(numedges(G), 3);
        end

    end

end

function d = dataDir()
%DATADIR   The package's data/ folder.
    d = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data');
end
