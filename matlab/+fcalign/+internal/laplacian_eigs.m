function [V, l, d] = laplacian_eigs(W, k, varargin)
% computes the eigendecomposition of the normalized Laplacian of a matrix
%
% Arguments
% ---------
% W : n x n numeric
%   affinity matrix
% k : 1 <= int <=n
%   rank of eigendecomposition to find
%
% Keyword Arguments
% -----------------
% center : logical
%   if true, centers laplacian before normalizing (default false)
% seed : o < int
%   random seed used for eigs (default 0)
%
% Returns
% -------
% V : k x n double
%   principal eigenvectors
% l : k x 1 double
%   principal eigenvalues
% d : n x 1 double
%   continuous degrees

parser = inputParser();
parser.addRequired('W', @(x) validateattributes(x, {'numeric'}, {'2d', 'square'}));
parser.addRequired('k', @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
parser.addParamValue('center', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
parser.addParamValue('seed', 0, @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
parser.parse(W, k, varargin{:});
inputs = parser.Results;

% further process inputs
[n, m] = size(W);

if k > n
    fcalign.Exception('InvalidAttribute', 'k must be less than or equal to the size of W');
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

% get eigen decomposition (assuming this is sorted in descending order)
[V, L] = eigs(Ms, k);
l = diag(L);
