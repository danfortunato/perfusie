% TEST_PERFUSIE   Run the perfusie unit tests.
%   Runs every Test*.m class in the tests/ folder and errors if any test
%   fails or any test file cannot be loaded, so that "mip test perfusie"
%   reports a failure. Tests that need the FMM are skipped (not failed)
%   when fmm2d is unavailable.

import matlab.unittest.TestSuite
import matlab.unittest.TestRunner

testdir = fullfile(fileparts(mfilename('fullpath')), 'tests');
files = dir(fullfile(testdir, 'Test*.m'));
suite = arrayfun(@(f) TestSuite.fromFile(fullfile(f.folder, f.name)), ...
    files, UniformOutput=false);
suite = [suite{:}];
if ( isempty(suite) )
    error('perfusie:noTests', 'No tests found in %s.', testdir);
end

runner = TestRunner.withTextOutput();
results = runner.run(suite);
disp(table(results))

nfailed = nnz([results.Failed]);
if ( nfailed > 0 )
    error('perfusie:testsFailed', '%d of %d tests failed.', nfailed, numel(results));
end
