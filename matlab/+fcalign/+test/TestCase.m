classdef (Abstract) TestCase < matlab.unittest.TestCase
    % test usage of a component of the fc align package

methods (Access = protected)

function data = create_data(this, varargin)
% creates an example data matrix
%
% creates an example data matrix, where each element
% is drawn from a uniform distribution
%
% Keyword Arguments
% -----------------
% n : 1 <= integer
%   number of nodes (default = 16)
% s : 1 <= integer
%   number of samples (default = 32)
% seed : 0 <= integer
%   seed for random number generator (default = 0)
%
% Returns
% -------
% data : n x s double
%   example data
    parser = inputParser();
    parser.addParamValue('n', 16, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.addParamValue('s', 32, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.parse(varargin{:});
    inputs = parser.Results;

    rng(inputs.seed);
    data = rand(inputs.n, inputs.s);
end

function matrix = create_correlation_matrix(this, varargin)
% creates a example correlation matrix
%
% creates an example correlation by creating a random example
% data matrix using this.create_data computing the pearson
% correlation coefficients and thresholding the values
%
% Keyword Arguments
% -----------------
% n : 1 <= integer
%   number of nodes (default = 16)
% s : 1 <= integer
%   number of samples (default = 32)
% threshold : -1 <= double <= 1
%   correlation threshold (default = 0.1)
% seed : 0 <= integer
%   seed for random number generator
%
% Returns
% -------
% data : n x s double
%   example correlation matrix
    parser = inputParser();
    parser.addParamValue('n', 16, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.addParamValue('s', 32, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.addParamValue('threshold', 0.1, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', -1, '<=', 1}));
    parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.parse(varargin{:});
    inputs = parser.Results;

    data = this.create_data('n', inputs.n, 's', inputs.s, 'seed', inputs.seed);
    matrix = fcalign.create_correlation_matrix(data, inputs.threshold);
end

function map = create_functional_map(this, varargin)
% creates a standard functional map
%
% creates an example functional map by creating an example
% correlation matrix with this.create_correlation_matrix
% and using default values to create the corresponding
% functional map
%
% Keyword Arguments
% -----------------
% n : 1 <= integer
%   number of nodes (default = 16)
% s : 1 <= integer
%   number of samples (default = 32)
% threshold : -1 <= double <= 1
%   correlation threshold (default = 0.1)
% seed : 0 <= integer
%   seed for random number generator (default = 0)
%
% Returns
% -------
% map : struct
%   example functional map
    parser = inputParser();
    parser.addParamValue('n', 16, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.addParamValue('s', 32, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.addParamValue('threshold', 0.1, @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', -1, '<=', 1}));
    parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.parse(varargin{:});
    inputs = parser.Results;

    matrix = this.create_correlation_matrix('n', inputs.n, 's', inputs.s, 'threshold', inputs.threshold, 'seed', inputs.seed);
    map = fcalign.create_functional_map(matrix);
end

function transform = create_transform(this, varargin)
% creates a p-dimension transform
%
% Keyword Arguments
% -----------------
% p : 1 <= integer
%   dimension of transform (default = 20)
% type : {'rotation', 'orthonormal', 'rigid', 'scale_orthonormal'}
%   type of transformation (default = 'rigid')
% seed : 0 <= integer
%   seed for random number generator (default = 0)
%
% Returns
% -------
% transform : p x p double
%   random transform
    parser = inputParser();
    parser.addParamValue('p', 20, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.addParamValue('type', 'rigid', @(x) validatestring(x, {'rotation', 'orthonormal', 'rigid', 'scale'}));
    parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'integer', 'scalar'}));
    parser.parse(varargin{:});
    inputs = parser.Results;

    rng(inputs.seed);

    % create a random orthonormal matrix
    orthonormal = orth(rand(p, p));

    % enforce determinant to be 1 to get rotation
    if strcmp(inputs.type, 'rotation') && strcmp(inputs.type, 'rigid')
        [V, L] = eig(rotation);
        if prod(diag(L)) < 0
            L(1,1) = -L(1, 1);
        end
        orthonormal = V*L*V';
    end

    % create a random translation
    translation = rand(p, 1);

    % create random transform
    transform = eye(p + 1, p + 1);
    switch inputs.type
        case 'rotation'
            transform(1:p, 1:p) = orthonormal;
        case 'orthonormal'
            transform(1:p, 1:p) = orthonormal;
        case 'rigid'
            transform(1:p, 1:p) = orthonormal;
            transform(1:p, p+1) = translation;
        case 'scale'
            scale = 0.1 + 2*rand;
            transform(1:p, 1:p) = scale*rotation;
            transform(1:p, p+1) = translation;
        otherwise
            throw(fcalign.Exception('TypeNotRecognized', ...
                sprintf('Transform type %s not recognized.', inputs.type)));
    end

end

function verify_equal_diffusion_map(this, actual, expected)
% verifies that the actual diffusion map is equal to the expected one
%
% Arguments
% ---------
% actual : struct
%   actual diffusion map
% expected : struct
%   expected diffusion map
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
% verifies that the actual functional map is equal to the expected one
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

end

end
