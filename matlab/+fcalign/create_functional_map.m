function map = create_functional_map(W, varargin)
% creates a functional map from the given correlation matrix
%
% *** This is a preliminary release, please don't distribute ***
%
% REQUIRED INPUTS
% W                 correlation matrix
%                   n x n double
%
% OPTIONAL INPUTS
% alpha             value defining member of alpha family
%                   0 <= double (default = 0)
% number_eigs       number of laplacian eigenvalues to compute
%                   1 <= integer <= n, t (default = 50)
% center_eigs       center laplacian
%                   logical (default = false)
% map_dimension     number of dimension in functional map
%                   1 <= integer <= n, t (default = 20)
% diffusion_time    diffusion time of functional map
%                   0 < double (default = 2)

% removed this for now
% delta             threshold of ratio between extreme eigenvalues
%                   double

% parse inputs
parser = inputParser();
parser.addRequired('W', @(x) validateattributes(x, {'double'}, {'2d'}));
parser.addParameter('alpha', [], @(x) validateattributes(x, {'double'}, {'nonnegative'}));
parser.addParameter('number_eigs', 50, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
parser.addParameter('center_eigs', false, @(x) validateattributes(x, {'logical'}));
parser.addParameter('map_dimension', 20, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
parser.addParameter('diffusion_time', 2, @(x) validateattributes(x, {'numeric'}, {'positive'}));
% parser.addParameter('delta', [], @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
parser.addParameter('verbose', false, @(x) validateattributes(x, {'logical'}));
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

[V, l, d] = fcalign.internal.laplacian_eigs(W, inputs.number_eigs, 'center', inputs.center_eigs);

if inputs.verbose
    fprintf('... done in %g s\n', toc(tic_id));
end

% create diffusion map
map = fcalign.create_diffusion_map(V, l, d, 'dim', inputs.map_dimension, 't', inputs.diffusion_time);
map.Ms_eig_vectors = V;
map.Ms_eig_values = l;