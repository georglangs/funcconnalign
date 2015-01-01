function matrix = create_correlation_matrix(data, threshold, varargin)
% creates a sparse correlation matrix from fMRI data
%
% Arguments
% ---------
% Y : n x t double
%   fMRI data for single subject (each row is a time series)
% threshold : -1 <= double <= 1
%   only keep correlations above this threshold
%
% Keyword Arguments
% -----------------
% keep_negatives : logical
%   if true, use absolute correlations (default = false)
% save_memory : logical
%   if true, create correlation matrix sparsely to save memory (default = true)
% min_density : 0 <= double <= 1
%   if negative, automatically choose better than chance (default = 0)
% verbose : logical
%   if true, prints updates as text (default = false)
%
% Returns
% -------
% matrix : struct
%   matrix and parameters, where
%   W: correlation matrix is an n' x n' double
%   non_zeros: non zero indices in data is an n x 1 integer
parser = inputParser();
parser.addRequired('data', @(x) validateattributes(x, {'numeric'}, {'2d'}));
parser.addRequired('threshold', @(x) validateattributes(x, {'double'}, {'scalar'}));
parser.addParamValue('save_memory', true, @(x) validateattributes(x, {'logical'}));
parser.addParamValue('min_density', 0, @(x) validateattributes(x, {'double'}));
parser.addParamValue('keep_negatives', false, @(x) validateattributes(x, {'logical'}));
parser.addParamValue('verbose', false, @(x) validateattributes(x, {'logical'}));
parser.parse(data, threshold, varargin{:});
inputs = parser.Results;

% store some inputs
matrix.threshold = inputs.threshold;
matrix.min_density = inputs.min_density;

% find zero rows
matrix.non_zeros = find(max(abs(data), [], 2) > 0);

if inputs.verbose
    fprintf('Found %u zero rows\n', numel(matrix.non_zeros));
end

% create sparse correlation coefficient matrix from non-zeros
if inputs.verbose
    fprintf('Generating correlation matrix ...\n');
    tic_id = tic;
end

data_non_zeros = data(matrix.non_zeros, :);
if inputs.save_memory
    matrix.W = fcalign.internal.sparse_corrcoef(data_non_zeros', inputs.threshold, 'keep_negatives', inputs.keep_negatives, 'verbose', inputs.verbose);
else
    matrix.W = corrcoef(data_non_zeros');
    to_keep = matrix.W > inputs.threshold;
    if inputs.keep_negatives
        to_keep = to_keep | matrix.W < inputs.threshold;
    end
    matrix.W = sparse(matrix.W .* to_keep);
    matrix.W(isnan(matrix.W)) = 0;
    clear W
end

if inputs.verbose
    fprintf('... done in %g s\n', toc(tic_id));
end

% calculate the degree of the nodes, and discard nodes with degree below threshold
if inputs.verbose
    fprintf('Sampling correlation matrix ...\n');
    tic_id = tic;
end

n_nodes = numel(matrix.non_zeros);

d_orig = full(sum(matrix.W, 2));
pi_orig = d_orig / sum(d_orig);

if inputs.min_density < 0
    nodes_kept = find(pi_orig > (1 / n_nodes));
else
    nodes_kept = find(pi_orig > inputs.min_density);
end

matrix.W = matrix.W(nodes_kept, nodes_kept);
matrix.node_indices = matrix.non_zeros(nodes_kept);

if inputs.verbose
    fprintf('... done in %g s\n', toc(tic_id));
end
