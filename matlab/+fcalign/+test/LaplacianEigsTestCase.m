classdef LaplacianEigsTestCase < fcalign.test.TestCase
    % Test usage of computing the laplacian eigenvalues.

methods (Test)

function test_default(obj)
    % Tests default usage of computing the laplacian eigenvalues.

    % define parameters
    n = 16;
    s = 29;
    k = 8;

    % make a positive definite matrix
    rng(0);
    data = rand(n, s);
    W = data * data';

    % create the expected eigenvalue decomposition
    expected_d = full(sum(W,2));
    inv_sqrt_D = spdiags(expected_d.^-0.5, 0, n, n);
    Ms = inv_sqrt_D * W * inv_sqrt_D;
    Ms = (Ms + Ms') / 2;
    rng(0);
    [expected_V, L] = eigs(Ms, k);
    expected_l = diag(L);

    % create correlation matrix using method
    [actual_V, actual_l, actual_d] = fcalign.laplacian_eigs(W, k);

    % verify that actual is equal to expected
    obj.verifyEqual(actual_V, expected_V, 'RelTol', 0.01);
    obj.verifyEqual(actual_l, expected_l, 'RelTol', 0.01);
    obj.verifyEqual(actual_d, expected_d, 'RelTol', 0.01);
end

end

end
