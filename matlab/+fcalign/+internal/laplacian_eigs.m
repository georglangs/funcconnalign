function [V, l, d] = laplacian_eigs(W, k, varargin)
% finds the eigendecomposition of the normalized Laplacian of a matrix
%
% [V, l, d] = laplacian_eigs(W, k, center)
%
% REQUIRED INPUTS
% W         affinity matrix
%           n x n double
% k         rank of eigendecomposition to find
%           1 <= integer <= n
%
% OPTIONAL INPUTS
% center    if true, centers laplacian before normalizing
%           logical (default false)
% seed      random seed used for eigs
%           0 < integer (default 0)
%
% OUTPUTS
% V         principal eigenvectors
%           k x n double
% l         principal eigenvalues
%           k x 1 double
% d         continuous degrees
%           n x 1 double

% parse inputs
parser = inputParser();
parser.addRequired('W', @(x) validateattributes(x, {'numeric'}, {'2d'}));
parser.addRequired('k', @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
parser.addParameter('center', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
parser.addParameter('seed', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
parser.parse(W, k, varargin{:});
inputs = parser.Results;

% further process inputs
[n, m] = size(W);

if n ~= m
    error('fcalign:laplacian_eigs:W_is_not_square', 'W must be a square matrix');
end

if k > n
    error('fcalign:laplacian_eigs:k_is_too_big', 'k must be less than or equal to the size of W');
end

% reseed random number generator
rng(inputs.seed);

% compute laplacian and compute sums
if inputs.center
    W = W - diag(diag(W));
end
d = full(sum(W,2));
inv_sqrt_D = spdiags(d.^-0.5, 0, n, n);

% compute normalized laplacian
Ms = inv_sqrt_D * W * inv_sqrt_D;

% ensure symmetry (this is a hack to get round numerical precision lost above)
Ms = (Ms + Ms') / 2;

% get eigen decomposition (this should be sorted in descending order)
[V, L] = eigs(Ms, k);
l = diag(L);