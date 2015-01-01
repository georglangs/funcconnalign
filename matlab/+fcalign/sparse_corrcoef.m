function [s, r, c] = sparse_corrcoef(X, threshold, varargin)
% creates a sparse correlation coefficient matrix from a data matrix and a threshold
%
% Arguments
% ---------
% X : n x p numeric
%   data matrix where each row is an observation
% threshold : numeric
%   lower threshold on correlation coefs or significance values
%
% Keyword Arguments
% -----------------
% threshold_type : {'correlation', 'significance'}
%   what to threshold, where
%   'correlation' corresponds to Pearson correlation coefficient
%   'significance' corresponds to Bonferroni corrected p values
%   (default = 'correlation')
% keep_negatives : boolean
%   threshold absolute values (default = false)
% nonzero_prop : 0 < double <= 1
%   maximum proportion of non-zero values we expect to get (default = 0.02)
% growth_rate : 0 < double
%   fraction of vector size to expand by when number of non-zeros exceeds (default = 0.25)
% verbose : boolean
%   if true, display messages about progress (default = true)
%
% Returns
% -------
% s : {n x n | nnz x 1} double
%   if only s is output, then this is the sparse correlation matrix, o.w. this is the dense vector of values exceeding the threshold
% r : nnz x 1 integer
%   row indices of nonzero matrix values
% c : nnz x 1 integer
%   column indices of nonzero matrix values

parser = inputParser;
parser.addRequired('X', @(x) validateattributes(x, {'numeric'}, {'2d'}));
parser.addRequired('threshold', @(x) validateattributes(x, {'numeric'}, {'scalar'}));
parser.addParamValue('threshold_type', 'correlation', @(x) validatestring(x, {'correlation', 'significance'}));
parser.addParamValue('keep_negatives', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
parser.addParamValue('nonzero_prop', 0.02, @(x) validateattributes(x, {'double'}, {'scalar', '>', 0, '<=', 1}));
parser.addParamValue('growth_rate', 0.25, @(x) validateattributes(x, {'double'}, {'scalar', '>', 0, '<=', 1}));
parser.addParamValue('verbose', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
parser.parse(X, threshold, varargin{:});
inputs = parser.Results;

% further process inputs
[n, p] = size(X);
inputs.nonzero_max = ceil(p^2 * inputs.nonzero_prop);
growth = ceil(inputs.growth_rate * inputs.nonzero_max);

% initialize index and value vectors
r = zeros(inputs.nonzero_max, 1);
c = zeros(inputs.nonzero_max, 1);
s = zeros(inputs.nonzero_max, 1);

% note that correlation coef matrix always has main diagonal of 1s (we use
% 0.5 because later we add the matrix in lower triangular form to its
% transpose to produce the actual matrix)
r(1:p) = 1:p;
c(1:p) = 1:p;
s(1:p) = 0.5*ones(p, 1);

% initially find all standard deviations and means
means = mean(X)';
stds = std(X)';

% zero mean all observations once
zmX = X - repmat(means', n, 1);

normalizer = 1 / (n-1);

% find correlations by variable, building lower triangular corr coef matrix
nnz = p;
if inputs.verbose
    fprintf('Computing sparse corrcoef: %03u%%', 0);
end
for i = 1:(p-1)

    % the unnnormalized covariance calc for the ith variable with the other n-i variables involves multiplying the matrix of observations of these other variables with the column of observations for the ith variable
    unnorm_covs = zmX(:,i+1:end)' * zmX(:,i);

    % normalize by constant to get covariance, then element wise by vector of standard deviation products to get correlation coefs
    cs = normalizer * unnorm_covs ./ (stds(i+1:end) * stds(i));

    % update progress
    if inputs.verbose
        fprintf('\b\b\b\b%03u%%', floor(i / p * 100));
    end

    % find the values that should be nonzero and their inds
    switch inputs.threshold_type
    case 'significance'
        n_dofs = n - 2;
        n_comparisons = p*(p-1)/2;

        ts = cs .* (n_dofs ./ (1 - cs.^2)).^0.5;

        if inputs.keep_negatives
            t_thresh = tinv(1 - 0.5*inputs.threshold / n_comparisons, n_dofs);
            nonzero_inds = find((ts > t_thresh) | (ts < -t_thresh));
        else
            t_thresh = tinv(1 - inputs.threshold / n_comparisons, n_dofs);
            nonzero_inds = find(ts > t_thresh);
        end
    case 'correlation'
        if inputs.keep_negatives
            nonzero_inds = find((cs > inputs.threshold) | (cs < -inputs.threshold));
        else
            nonzero_inds = find(cs > inputs.threshold);
        end
    otherwise
        throw(fcalign.Exception('UnrecognizedValue', 'Threshold type %s is not recognized', inputs.threshold_type));
    end

    nonzero_cs = cs(nonzero_inds);

    % if adding these takes us over max, then add some more space
    % (proportion to original max)
    nnz_cs = length(nonzero_cs);
    if (nnz + nnz_cs) > length(s)
        r = [r; zeros(growth, 1)];
        c = [c; zeros(growth, 1)];
        s = [s; zeros(growth, 1)];
    end

    % inds in dense vector
    inds = (nnz+1):(nnz+nnz_cs);
    % row inds in matrix (note we're building lower triangular, so shift these by i)
    r(inds) = i + nonzero_inds;
    % this is correlation for ith variable
    c(inds) = i;

    s(inds) = nonzero_cs;

    nnz = nnz + nnz_cs;
end

if (nargout > 1)
    r(s == 0) = [];
    c(s == 0) = [];
    s(s == 0) = [];
else
    upper_tri = sparse(r(1:nnz), c(1:nnz), s(1:nnz), p, p, 2*nnz - p);
    s = upper_tri + upper_tri';
end

if inputs.verbose
    fprintf('\b\b\b\b%03u%%', 100);
    fprintf('\n');
end
