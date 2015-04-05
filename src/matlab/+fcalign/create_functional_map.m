function map = create_functional_map(W, varargin)
% Creates a functional map from the given correlation matrix.
%
% Arguments
% ---------
% W : n x n double
%   correlation matrix
%
% Keyword Arguments
% -----------------
% alpha : 0 <= double
%   value defining member of alpha family (default = 0)
% number_eigs : 1 <= integer <= n, t
%   number of laplacian eigenvalues to compute (default = 50)
% center_eigs : logical
%   center laplacian matrix before taking eigendecomposition (default = false)
% map_dimension : 1 <= integer <= n, t
%   number of dimension in functional map (default = 20)
% diffusion_time : 0 < double
%   diffusion time of functional map (default = 2)

% removed this for now
% delta             threshold of ratio between extreme eigenvalues
%                   double
parser = inputParser();
parser.addRequired('W', @(x) validateattributes(x, {'double'}, {'2d'}));
parser.addParamValue('alpha', [], @(x) validateattributes(x, {'double'}, {'nonnegative'}));
parser.addParamValue('number_eigs', 50, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
parser.addParamValue('center_eigs', false, @(x) validateattributes(x, {'logical'}));
parser.addParamValue('map_dimension', 20, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
parser.addParamValue('diffusion_time', 2, @(x) validateattributes(x, {'numeric'}, {'positive'}));
% parser.addParamValue('delta', [], @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
parser.addParamValue('verbose', false, @(x) validateattributes(x, {'logical'}));
parser.parse(W, varargin{:});
inputs = parser.Results;

% process inputs further
n = size(W, 1);

if inputs.number_eigs > n
    error('fcalign:create_functional_map:number_eigs_too_large', 'Number of eigenvalues must not be greater than size of correlation matrix');
end

if inputs.map_dimension > inputs.number_eigs
    error('fcalign:create_functional_map:map_dimension_too_large', 'Map dimension must not be greater than number of eigenvalues');
end

if inputs.verbose
    fprintf('... done\n');
end

% normalize correlation matrix by power of density according to family parameter
if inputs.alpha > 0
    if inputs.verbose
        fprintf('Computing member of alpha family ...\n');
        tic_id = tic;
    end

    d = full(sum(W, 2));
    inv_D_alpha = spdiags(d.^-inputs.alpha, 0, n, n);
    W = inv_D_alpha*W*inv_D_alpha;

    if inputs.verbose
        fprintf('... done in %g s\n', toc(tic_id));
    end
end

% get eigenvectors of normalized laplacian
if inputs.verbose
    fprintf('Taking eigendecomposition of normalized Laplacian ...\n');
    tic_id = tic;
end

[V, l, d] = fcalign.laplacian_eigs(W, inputs.number_eigs, 'center', inputs.center_eigs);

if inputs.verbose
    fprintf('... done in %g s\n', toc(tic_id));
end

% create diffusion map
map = fcalign.create_diffusion_map(V, l, d, 'dim', inputs.map_dimension, 't', inputs.diffusion_time);
map.Ms_eig_vectors = V;
map.Ms_eig_values = l;
