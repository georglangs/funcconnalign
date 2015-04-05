classdef CreateCorrelationMatrixTestCase < fcalign.test.TestCase
   % Test usage of creating a sparse correlation matrix.

methods (Test)

function test_default(this)
    % Tests default usage of creating a correlation matrix.

    % define parameters
    n = 16;
    s = 32;
    zero_inds = [3, 8];
    threshold = 0.2;

    % define some small random fmri data
    data = this.create_data(n, s);

    % make some of the rows be zeros
    data(zero_inds, :) = 0;

    % create the expected sparse correlation matrix
    expected.non_zeros = setdiff(1:n, zero_inds)';
    non_zero_data = data(expected.non_zeros, :);
    expected.W = corrcoef(non_zero_data');
    expected.W = expected.W .* (expected.W > threshold);
    expected.W = sparse(expected.W);

    % create correlation matrix using method
    actual = fcalign.create_correlation_matrix(data, threshold);

    % verify that expected correlation
    this.verifyEqual(actual.W, expected.W, 'RelTol', 0.01);
    this.verifyEqual(actual.non_zeros, expected.non_zeros);
end

end

end
