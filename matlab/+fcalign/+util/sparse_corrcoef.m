function [s, r, c] = sparse_corrcoef(X, threshold, varargin)
% creates a sparse correlation coefficient matrix from a data matrix and a threshold
%
% function [s, r, c] = sparse_corrcoef(X, threshold, varargin)
%
% REQUIRED INPUTS
% X                 data matrix where each row is an observation
%                   n x p double
% threshold         lower threshold on correlation coefs or significance values
%                   double
%
% OPTIONAL INPUTS
% threshold_type    what to threshold
%                   string from:
%                       'correlation': Pearson correlation coefficient (default)
%                       'significance': Bonferroni corrected p-value of Pearson correlation coef
% keep_negatives    threshold absolute values
%                   boolean (default = false)
% nonzero_max       maximum number of non-zero values we expect to get
%                   integer (default = ???)
% growth_rate       fraction of vector size to expand by when number of non-zeros exceeds
%                   0 < double
% verbose           if true, display messages about progress
%                   boolean (true)
%
% OUTPUTS
% s                 if only s is output, then this is the sparse correlation matrix
%                   sparse nxn double
%                   o.w. this is the dense vector of values exceeding the threshold
%                   nnz x 1 double
% r                 row indices of nonzero matrix values
%                   nnz x 1 integer
% c                 column indices of nonzero matrix values
%                   nnz x 1 integer

[n, p] = size(X);

parser = inputParser;
parser.addRequired('X', @(x) ismatrix(x) && isnumeric(x));
parser.addRequired('threshold', @(x) isscalar(x) && isnumeric(x));
parser.addParameter('threshold_type', 'correlation', @(x) validatestring(x, {'correlation', 'significance'}));
parser.addParameter('keep_negatives', false, @(x) isscalar(x) && islogical(x));
parser.addParameter('nonzero_max', , @(x) isscalar(x) && isnumeric(x) && ~mod(x, 1));
parser.addParameter('growth_rate', 0.25, @(x) isscalar(x) && isnumeric(x) && (x > 0));
parser.addParameter('verbose', false, @(x) isscalar(x) && islogical(x));
parser.parse(X, threshold, varargin{:});
inputs = parser.Results;

% if we're not given the expected max number of nonzeros, then assume sparsity won't exceed 1/50
if ~exist('nonzero_max', 'var') || isempty(nonzero_max)
    nonzero_max = ceil((p^2)/50);
    if verbose
        disp(['using auto nonzero_max ', num2str(nonzero_max)]);
    end
end

r = zeros(nonzero_max, 1);
c = zeros(nonzero_max, 1);
s = zeros(nonzero_max, 1);

% we grow these these vectors if they're full at a rate proportional to the original max
growth_rate = 0.5;
growth = int32(growth_rate*nonzero_max);

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

normalizer = 1/(n-1);

% find correlations by variable, building lower triangular corr coef matrix
nnz = p;
for i=1:(p-1)

    if verbose
        if (~mod(i,500))
            disp(['calculated for ', num2str(i), ' variables']);
        end
    end

    % the unnnormalized covariance calc for the ith variable with the other n-i variables involves multiplying the matrix of observations of these other variables with the column of observations for the ith variable
    unnorm_covs = zmX(:,i+1:end)'*zmX(:,i);

    % normalize by constant to get covariance, then element wise by vector of standard deviation products to get correlation coefs
    cs = normalizer * unnorm_covs ./ (stds(i+1:end)*stds(i));

    % find the values that should be nonzero and their inds
    if (strcmp(thresh_type, 'significance'))

        n_dofs = n - 2;
        n_comparisons = p*(p-1)/2;

        ts = cs .* (n_dofs ./ (1 - cs.^2)).^0.5;

        if (keep_negatives)
            t_thresh = tinv(1 - 0.5*thresh/n_comparisons, n_dofs);
            nonzero_inds = find((ts > t_thresh) | (ts < -t_thresh));
        else
            t_thresh = tinv(1 - thresh/n_comparisons, n_dofs);
            nonzero_inds = find(ts > t_thresh);
        end

    else
        if (keep_negatives)
            nonzero_inds = find((cs > thresh) | (cs < -thresh));
        else
            nonzero_inds = find(cs > thresh);
        end
    end

    nonzero_cs = cs(nonzero_inds);

    % if adding these takes us over max, then add some more space
    % (proportion to original max)
    nnz_cs = length(nonzero_cs);
    if ( nnz + nnz_cs > length(s) )
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
    % knock off zero elements
    r(s == 0) = [];
    c(s == 0) = [];
    s(s == 0) = [];

else
    upper_tri = sparse(r(1:nnz), c(1:nnz), s(1:nnz), p, p, 2*nnz - p);
    s = upper_tri + upper_tri';
end
