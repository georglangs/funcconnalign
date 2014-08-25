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
% delta             threshold of ratio between extreme eigenvalues
%                   double
%
% Georg Langs
% 23.9.2011
%
% Edited by Andy Sweet on 23.8.2014

%% parse inputs
if p.verbose
    fprintf('Parsing inputs...\n');
end

parser = inputParser();
parser.addRequired('W', @(x) validateattributes(x, {'double'}, {'2d'}));
parser.addParameter('alpha', [], @(x) validateattributes(x, {'double'}, {'nonnegative'}));
parser.addParameter('number_eigs', 50, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
parser.addParameter('center_eigs', false, @(x) validateattributes(x, {'logical'}));
parser.addParameter('map_dimension', 20, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
parser.addParameter('diffusion_time', 2, @(x) validateattributes(x, {'numeric'}, {'positive'}));
parser.addParameter('delta', [], @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
parser.addParameter('verbose', true, @(x) validateattributes(x, {'logical'}));
parser.parse(data, varargin{:});
p = parser.Results;

% validate and process inputs further
n = size(W, 1);

if p.number_eigs > n
    error('fcalign:create_functional_map:number_eigs_too_large', 'Number of eigenvalues must not be greater than size of correlation matrix');
end

if p.map_dimension > p.number_eigs
    error('fcalign:create_functional_map:map_dimension_too_large', 'Map dimension must not be greater than number of eigenvalues');
end

if p.verbose
    fprintf('... done\n');
end

% store in map
map.p = p;

%% compute degrees of correlation matrix
map.d = full(sum(W, 2));

%% normalize correlation matrix by power of density according to family parameter
if p.alpha > 0
    if p.verbose
        fprintf('Computing member of alpha family ...\n');
        ticId = tic;
    end

    inv_D_alpha = spdiags(map.d.^-p.alpha, 0, n, n);
    W = inv_D_alpha*W*inv_D_alpha;

    if p.verbose
        fprintf('... done in %g s\n', toc(ticId));
    end
end

%% get eigenvectors of normalized laplacian
if p.verbose
    fprintf('Taking eigendecomposition of normalized Laplacian ...\n');
    ticId = tic;
end

[map.Ms_eig_vectors, map.Ms_eig_values] = laplacian_eigs(W, p.number_eigs, p.center_eigs);

if p.verbose
    fprintf('... done in %g s\n', toc(ticId));
end

%% create functional map
[map.Gamma, determined_param] = create_diffusion_map(map.Ms_eig_vectors, map.Ms_eig_values, map.d, p.map_dimension, p.diffusion_time, p.delta);

if isempty(p.diffusion_time)
    if verbose
        fprintf('Determined diffusion time = %g', determined_param);
    end
    map.DeterminedDiffusionTime = determined_param;
elseif isempty(p.map_dimension)
    if verbose
        fprintf('Determined map dimension = %u', determined_param);
    end
    map.DeterminedMapDimension = determined_param;
elseif isempty(p.delta)
    if verbose
        fprintf('Determined delta = %g', determined_param);
    end
    map.DeterminedDelta = determined_param;
end

map.Delta = p.delta;
map.DiffusionTime = p.diffusion_time;
map.MapDimension = p.map_dimension;
