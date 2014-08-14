function [s,r,c] = atlas_sparse_cross_corrcoef(X, Y, thresh, thresh_type, keep_negatives, nzmax, growth_rate, verbose)
% ATLAS_SPARSE_CROSS_CORRCOEF creates a sparse correlation coefficient matrix given two data matrices X and Y and a lower bound thresh
% 
% INPUTS
% X                 n x p: data matrix where each row is an observation (required)
% Y                 n x p: data matrix where each row is an observation (required)
% thresh            lower threshold on correlation coefs or significance values (required)
% thresh_type       what to threshold
%                       'correlation': Pearson correlation coefficient (default)
%                       'significance': Bonferroni corrected p-value of Pearson correlation coef
% keep_negatives    flag: threshold absolute values 
% nzmax             the maximum number of non-zero values we expect to get
% growth_rate       fraction of current vector size to expand by when nnzs exceeds it
% verbose           flag: display messages about progress
%
% OUTPUTS
% s                 n x n: if only s is output, then this is the sparse correlation matrix
%                   nnz x 1: o.w. this is the dense vector of values exceeding the threshold
% r                 nnz x 1: row inds of nonzero matrix values
% c                 nnz x 1: column inds of nonzero matrix values
%
% Andy Sweet 02/12/2010

[n_x, p_x] = size(X);
[n_y, p_y] = size(Y);

if ( n_x ~= n_y )
    error(['Number of observations in X, ', num2str(n_x), ', does not equal number of observations in Y, ', num2str(n_y), '.']);
end
n = n_x;
p = (p_x + p_y)/2;

if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end

if ~exist('growth_rate', 'var') || isempty(growth_rate)
    growth_rate = 0.25;
end

if ~exist('keep_negatives', 'var') || isempty(keep_negatives)
    keep_negatives = false;
end

if ~exist('thresh_type', 'var') || isempty(thresh_type)
    thresh_type = 'distance';
end

% if we're not given the expected max number of nonzeros, then assume sparsity won't exceed 1/50
if ~exist('nzmax', 'var') || isempty(nzmax)
    nzmax = ceil((p^2)/50);
    if verbose
        disp(['using auto nzmax ', num2str(nzmax)]);
    end
end

r = zeros(nzmax, 1);
c = zeros(nzmax, 1);
s = zeros(nzmax, 1);

% we grow these these vectors if they're full at a rate proportional to the original max
growth_rate = 0.25;
growth = int32(growth_rate*nzmax);

% initially find all standard deviations and means
means_x = mean(X)';
stds_x = std(X)';
means_y = mean(Y)';
stds_y = std(Y)';

% zero mean all observations once
zmX = X - repmat(means_x', n_x, 1);
zmY = Y - repmat(means_y', n_y, 1);

normalizer = 1/(n-1);

% find correlations by variable, building lower triangular corr coef matrix
nnz = 0;
for i=1:p_x
    
    % unnormalized cross covariances between ith x variable and all ys
    unnorm_covs = zmX(:,i)'*zmY;
    
    % normalize by constant to get covariance, then element wise by vector of standard deviation products to get correlation coefs
    cs = normalizer * unnorm_covs ./ (stds_x(i)*stds_y');
    
    % find the values that should be nonzero and their inds
    if (strcmp(thresh_type, 'significance'))
        
        n_dofs = n - 2;
        n_comparisons = p_x*p_y;
        
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
    r(inds) = i; 
    % this is correlation for ith variable
    c(inds) = nonzero_inds;
    s(inds) = nonzero_cs;
    
    nnz = nnz + nnz_cs;
end

if (nargout > 1)
    % knock off zero elements
    r(s == 0) = [];
    c(s == 0) = [];
    s(s == 0) = [];
else
    s = sparse(r(1:nnz), c(1:nnz), s(1:nnz), p_x, p_y, nnz);
end
