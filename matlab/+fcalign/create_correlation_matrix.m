function matrix = create_correlation_matrix(data, varargin)
% creates a sparse correlation matrix from fMRI data
%
% *** This is a preliminary release, please don't distribute ***
%
% REQUIRED INPUTS
% Y                 fMRI data for single subject (each row is a time series)
%                   n x t double
%
% OPTIONAL INPUTS
% threshold         only keep correlations above this threshold
%                   double
% keep_negatives    if true, use absolute correlations
%                   logical (default = false)
% save_memory       if true, create correlation matrix sparsely to save memory
%                   logical (default = true)
% min_density       if negative, automatically choose better than chance
%                   0 <= double <=1 (default = -1)
%
% OUTPUTS
% matrix            matrix and parameters
%                   struct:
%                       W: correlation matrix
%                       n' x n' double
%                       non_zeros: non zero indices in data
%                       n x 1 integer
%
%
% Georg Langs
% 23.9.2011
%
% Edited by Andy Sweet on 23.8.2014

%% parse inputs
parser = inputParser();
parser.addRequired('data', @(x) validateattributes(x, {'struct'}));
parser.addParameter('save_memory', true, @(x) validateattributes(x, {'logical'}));
parser.addParameter('min_density', -1, @(x) validateattributes(x, {'double'}));
parser.addParameter('threshold', 0.4, @(x) validateattributes(x, {'double'}));
parser.addParameter('keep_negatives', false, @(x) validateattributes(x, {'logical'}));
parser.addParameter('verbose', true, @(x) validateattributes(x, {'logical'}));
parser.parse(data, varargin{:});
p = parser.Results;

%% store some inputs
matrix.threshold = p.threshold;
matrix.min_density = p.min_density;

%% find zero rows
matrix.non_zeros = find(max(abs(data), [], 2) > 0);

if p.verbose
    fprintf('Found %u zero rows', numel(matrix.non_zeros));
end

%% create sparse correlation coefficient matrix from non-zeros
if verbose
    fprintf('Generating correlation matrix ...\n');
    ticId = tic;
end

if p.save_memory
    matrix.W = sparse_corrcoef(data(data.non_zeros, :)', p.threshold, 'correlation', p.keep_negatives, p.verbose);
else
    W = corrcoef(data.BOLD(data.non_zeros, :)');
    matrix.W = sparse(W.*(W > p.threshold));
    clear W
    matrix.W(isnan(matrix.W)) = 0;
end

if p.verbose
    fprintf('... done in %g s\n', toc(ticId));
end

%% calculate the degree of the nodes, and discard nodes with degree below threshold
if verbose
    fprintf('Sampling correlation matrix ...\n');
    ticId = tic;
end

n_nodes = numel(matrix.non_zeros);

d_orig = full(sum(matrix.W, 2));
pi_orig = d_orig / sum(d_orig);

if p.min_density < 0
    nodes_kept = find(pi_orig > (1 / n_nodes));
else
    nodes_kept = find(pi_orig > p.min_density);
end

matrix.W = matrix.W(nodes_kept, nodes_kept);
matrix.node_indices = data.non_zeros(nodes_kept);

if p.verbose
    fprintf('... done in %g s', toc);
end
