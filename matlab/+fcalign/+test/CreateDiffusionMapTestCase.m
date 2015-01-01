classdef CreateDiffusionMapTestCase < fcalign.test.TestCase
    % tests different usages of creating a diffusion map

methods (Test)

function test_default(this)
% tests default usage of creating a functional map

    % define parameters
    n = 16;
    s = 32;
    k = 7;
    p = 5;
    t = 1;

    % create a correlation matrix
    matrix = this.create_correlation_matrix('n', n, 's', s);

    % get the laplacian eigs
    [V, l, d] = fcalign.laplacian_eigs(matrix.W, k);

    % create the expected sparse correlation matrix
    norm_factor = sqrt(sum(d));
    inv_sqrt_D = spdiags(d.^-0.5, 0, n, n);
    sqrt_D = spdiags(d.^0.5, 0, n, n);
    expected.Phi = sqrt_D*V(:, 1:p) / norm_factor;
    expected.Psi = norm_factor * inv_sqrt_D*V(:, 1:p);

    if expected.Psi(1, 1) < 0
    expected.Phi = expected.Phi * -1;
    expected.Psi = expected.Psi * -1;
    end

    Lt = spdiags(l(1:p).^t, 0, p, p);
    expected.Gamma = expected.Psi*Lt;

    expected.diffusion_time = t;
    expected.map_dimension = p;
    expected.delta = (l(p) / l(2))^t;
    expected.determined_param = 'delta';

    % create correlation matrix using method
    actual = fcalign.create_diffusion_map(V, l, d, 'dim', p, 't', t);

    % verify that expected correlation
    this.verify_equal_diffusion_map(actual, expected);
end

end

end
