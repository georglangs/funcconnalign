%% tests the sorted eigs function
function tests = sorted_eigs_test
    tests = functiontests(localfunctions);
end

%% test functions
function testDefault(testCase)

    % create normalized eigenvectors and eigenvalues
    v1 = [3; 9; -1];
    v2 = [-2.4; 3.2; 0];
    v3 = [0; -0.2; 1];
    v4 = [-5; 0; -5];
    V = [v1, v2, v3, v4];
    V = bsxfun(@rdivide, V, sqrt(sum(V.^2, 2)));
    l = [5, -1, 0.2];
    L = diag(l);

    % create matrix from these
    A = V*L*V';

    % compute eigendecomposition
    [actual.V, actual.l] = fcalign.util.sorted_eigs(A);

    % verify that it's equal within some tolerance
    testCase.verifyEqual(actual.V, , 'RelTol', 0.01);
    testCase.verifyEqual(actual.V, L, 'RelTol', 0.01);

end
