classdef SparseCorrcoefTestCase < fcalign.test.TestCase
    % Test usage of creating a sparse correlation matrix.

methods (Test)

function test_default(obj)
    % Tests default usage of the sparse corrcoef function.

    % define parameters
    n = 16;
    p = 29;
    threshold = 0.2;

    % reseed rng
    rng(0);

    % define some small random data
    data = rand(n, p);

    % create the expected sparse correlation matrix
    expected = corrcoef(data');
    expected = expected .* (expected > threshold);
    expected = sparse(expected);

    % create correlation matrix using method
    actual = fcalign.sparse_corrcoef(data', threshold);

    % verify that actual is equal to expected
    obj.verifyEqual(actual, expected, 'RelTol', 0.01);
end

end

end
