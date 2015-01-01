classdef CreateFunctionalMapTestCase < fcalign.test.TestCase
    % tests different usages of creating a functional map

methods (Test)

function test_default(this)
% tests default usage of creating a functional map

    % define parameters
    n = 16;
    s = 32;
    k = 9;
    dim = 6;
    t = 2;

    % create a correlation matrix
    matrix = this.create_correlation_matrix('n', n, 's', s);

    % create the expected map
    [V, l, d] = fcalign.laplacian_eigs(matrix.W, k);
    expected = fcalign.create_diffusion_map(V, l, d, 'dim', dim, 't', t);
    expected.Ms_eig_vectors = V;
    expected.Ms_eig_values = l;

    % create the actual map
    actual = fcalign.create_functional_map(matrix.W, ...
            'number_eigs', k, 'map_dimension', dim, ...
            'diffusion_time', t);

    % verify that actual functional map is equal to expected one
    this.verify_equal_functional_map(actual, expected);
end

end

end
