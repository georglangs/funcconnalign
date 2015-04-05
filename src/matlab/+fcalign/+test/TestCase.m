classdef (Abstract) TestCase < matlab.unittest.TestCase
    % Test usage of a component of the fcalign package.

methods (Access = protected)

function data = create_data(this, n, s, varargin)
    % Creates an example data matrix.
    %
    % Creates an example data matrix, where each element
    % is drawn from a uniform distribution.
    %
    % Arguments
    % ---------
    % n : 1 <= integer
    %   Dimension of data.
    % s : 1 <= integer
    %   Number of samples (default = 32).
    %
    % Keyword Arguments
    % -----------------
    % seed : 0 <= integer
    %   Seed for random number generator (default = 0).
    %
    % Returns
    % -------
    % data : n x s double
    %   Example data.

    parser = inputParser();
    parser.addRequired('n', @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar', 'positive'}));
    parser.addRequired('s', @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar', 'positive'}));
    parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar', 'nonnegative'}));
    parser.parse(n, s, varargin{:});
    inputs = parser.Results;

    rng(inputs.seed);
    data = rand(n, s);
end

function matrix = create_correlation_matrix(this, n, s, varargin)
    % Creates a example correlation matrix
    %
    % Creates an example correlation by creating a random example
    % data matrix using this.create_data computing the pearson
    % correlation coefficients and thresholding the values.
    %
    % Arguments
    % ---------
    % n : 1 <= integer
    %   Dimension of data.
    % s : 1 <= integer
    %   Number of samples.
    %
    % Keyword Arguments
    % -----------------
    % threshold : -1 <= double <= 1
    %   Correlation threshold (default = 0.1).
    % seed : 0 <= integer
    %   Seed for random number generator (default = 0).
    %
    % Returns
    % -------
    % data : n x s double
    %   Example correlation matrix.

    parser = inputParser();
    parser.addRequired('n', @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar', 'positive'}));
    parser.addRequired('s', @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar', 'positive'}));
    parser.addParamValue('threshold', 0.1, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', -1, '<=', 1}));
    parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar', 'nonnegative'}));
    parser.parse(n, s, varargin{:});
    inputs = parser.Results;

    data = this.create_data(n, s, 'seed', inputs.seed);
    matrix = fcalign.create_correlation_matrix(data, inputs.threshold);
end

function map = create_functional_map(this, n, s, varargin)
    % Creates a random functional map.
    %
    % Creates a random functional map by creating an random
    % correlation matrix with this.create_correlation_matrix
    % and using default values to create the corresponding
    % functional map.
    %
    % Keyword Arguments
    % -----------------
    % n : 1 <= integer
    %   Dimension of data.
    % s : 1 <= integer
    %   Number of samples.
    % threshold : -1 <= double <= 1
    %   Correlation threshold (default = 0.1).
    % seed : 0 <= integer
    %   Seed for random number generator (default = 0).
    %
    % Returns
    % -------
    % map : struct
    %   Random functional map.

    parser = inputParser();
    parser.addRequired('n', @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.addRequired('s', @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.addParamValue('threshold', 0.1, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', -1, '<=', 1}));
    parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar', 'nonnegative'}));
    parser.parse(n, s, varargin{:});
    inputs = parser.Results;

    matrix = this.create_correlation_matrix(n, s, 'threshold', inputs.threshold, 'seed', inputs.seed);
    map = fcalign.create_functional_map(matrix);
end

function transform = create_transform(this, n, varargin)
    % Creates a random n-dimensional transform.
    %
    % Arguments
    % ---------
    % n : 1 <= integer
    %   Dimension of transform.
    %
    % Keyword Arguments
    % -----------------
    % type : {'rotation', 'orthonormal', 'rigid', 'scale_orthonormal'}
    %   Type of transform (default = 'rigid').
    % seed : 0 <= integer
    %   Seed for random number generator (default = 0).
    %
    % Returns
    % -------
    % transform : n x n double
    %   Random transform.

    parser = inputParser();
    parser.addRequired('n', @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.addParamValue('type', 'rigid', @(x) validatestring(x, {'rotation', 'orthonormal', 'rigid', 'scale'}));
    parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar', 'nonnegative'}));
    parser.parse(n, varargin{:});
    inputs = parser.Results;

    rng(inputs.seed);

    % create a random orthonormal matrix
    orthonormal = orth(rand(n, n));

    % enforce determinant to be 1 to get rotation
    if strcmp(inputs.type, 'rotation') && strcmp(inputs.type, 'rigid')
        [V, L] = eig(rotation);
        if prod(diag(L)) < 0
            L(1, 1) = -L(1, 1);
        end
        orthonormal = V * L * V';
    end

    % create a random translation
    translation = rand(n, 1);

    % create random transform
    transform = eye(n + 1, n + 1);
    switch inputs.type
        case 'rotation'
            transform(1:n, 1:n) = orthonormal;
        case 'orthonormal'
            transform(1:n, 1:n) = orthonormal;
        case 'rigid'
            transform(1:n, 1:n) = orthonormal;
            transform(1:n, n+1) = translation;
        case 'scale'
            scale = 0.1 + 2*rand;
            transform(1:n, 1:n) = scale*rotation;
            transform(1:n, n+1) = translation;
        otherwise
            throw(fcalign.Exception('TypeNotRecognized', ...
                'Transform type %s not recognized.', inputs.type));
    end

end

function verify_equal_diffusion_map(this, actual, expected)
    % Verifies that the actual diffusion map is equal to the expected one.
    %
    % Arguments
    % ---------
    % actual : struct
    %   Actual diffusion map.
    % expected : struct
    %   Expected diffusion map.

    this.assertClass(actual, 'struct');
    this.assertClass(expected, 'struct');

    this.verifyEqual(actual.Phi, expected.Phi, 'RelTol', 0.01);
    this.verifyEqual(actual.Psi, expected.Psi, 'RelTol', 0.01);
    this.verifyEqual(actual.Gamma, expected.Gamma, 'RelTol', 0.01);
    this.verifyEqual(actual.diffusion_time, expected.diffusion_time, 'RelTol', 0.01);
    this.verifyEqual(actual.map_dimension, expected.map_dimension, 'RelTol', 0.01);
    this.verifyEqual(actual.delta, expected.delta, 'RelTol', 0.01);
    this.verifyEqual(actual.determined_param, expected.determined_param, 'RelTol', 0.01);
end

function verify_equal_functional_map(this, actual, expected)
    % Verifies that the actual functional map is equal to the expected one.
    %
    % Arguments
    % ---------
    % actual : struct
    %   actual functional map
    % expected : struct
    %   expected functional map

    this.assertClass(actual, 'struct');
    this.assertClass(expected, 'struct');

    this.verifyEqual(actual.Ms_eig_vectors, expected.Ms_eig_vectors, 'RelTol', 0.01);
    this.verifyEqual(actual.Ms_eig_values, expected.Ms_eig_values, 'RelTol', 0.01);
    this.verify_equal_diffusion_map(actual, expected);
end

function verify_equal_transform(this, actual, expected)
    % Verifies that the actual transform is equal to the expected one.
    %
    % Arguments
    % ---------
    % actual : nxn numeric
    %   actual transform
    % expected : nxn numeric
    %   expected transform
   
    this.verifyEqual(actual, expected, 'RelTol', 0.01);
end

end

end
