function [s,r,c] = sparse_corrcoef_fast(X, thresh, thresh_type, keep_negatives, nzmax)
% sparse_corrcoef creates a sparse correlation coefficient matrix given a data matrix X and a coefficient lower bound
% As opposed to sparse_corrcoef this stores the values and their positions in the matrix densely, then creates the sparse matrix at the end
% 
% INPUTS
% X: data matrix where each row is an observation
% thresh: lower bound of correlation coefficients to keep (to achieve sparsity)
% nzmax: the maximum number of non-zero values we expect to get (optional)
%
% OUTPUTS
% coef_mat: correlation coefficient matrix
%

[n_obs, n_vars] = size(X);

if (nargin > 4)
    disp(['using nzmax ', num2str(nzmax)]);
else
    % if we're not given the expected max number of nonzeros, then assume
    % this has ratio 1/50
    nzmax = ceil((n_vars^2)/50);
    disp(['using auto nzmax ', num2str(nzmax)]);
end

r = zeros(nzmax, 1);
c = zeros(nzmax, 1);
s = zeros(nzmax, 1);

% we grow these these vectors if they're full at a rate proportional to the original max
growth_rate = 0.5;
growth = int32(growth_rate*nzmax);

% note that correlation coef matrix always has main diagonal of 1s (we use
% 0.5 because later we add the matrix in lower triangular form to its
% transpose to produce the actual matrix)
r(1:n_vars) = 1:n_vars;
c(1:n_vars) = 1:n_vars;
s(1:n_vars) = 0.5*ones(n_vars, 1);

% initially find all standard deviations and means
means = mean(X)';
stds = std(X)';

% zero mean all observations once
zmX = X - repmat(means', n_obs, 1);

normalizer = 1/(n_obs-1);

% find correlations by variable, building lower triangular corr coef matrix
nnz = n_vars;
for i=1:(n_vars-1)
    
    if (~mod(i,500))
        disp(['calculated for ', num2str(i), ' variables']);
    end
    
    % the unnnormalized covariance calc for the ith variable with the other n-i variables involves multiplying the matrix of observations of these other variables with the column of observations for the ith variable
    unnorm_covs = zmX(:,i+1:end)'*zmX(:,i);
    
    % normalize by constant to get covariance, then element wise by vector of standard deviation products to get correlation coefs
    cs = normalizer * unnorm_covs ./ (stds(i+1:end)*stds(i));
    
    % find the values that should be nonzero and their indices
    if (strcmp(thresh_type, 'significance'))
        
        n_dofs = n_obs - 2;
        n_comparisons = n_vars*(n_vars-1)/2;
        
        ts = cs .* (n_dofs ./ (1 - cs.^2)).^0.5;
        
        if (keep_negatives)
            t_thresh = tinv(1 - 0.5*thresh/n_comparisons, n_dofs);
            nonzero_indices = find((ts > t_thresh) | (ts < -t_thresh));
        else
            t_thresh = tinv(1 - thresh/n_comparisons, n_dofs);
            nonzero_indices = find(ts > t_thresh);
        end

    else
        
        if (keep_negatives)
            nonzero_indices = find((cs > thresh) | (cs < -thresh));
        else
            nonzero_indices = find(cs > thresh);
        end
        
    end
    
    
    nonzero_cs = cs(nonzero_indices);
    
    % if adding these takes us over max, then add some more space
    % (proportion to original max)
    nnz_cs = length(nonzero_cs);
    if ( nnz + nnz_cs > length(s) )
        r = [r; zeros(growth, 1)];
        c = [c; zeros(growth, 1)];
        s = [s; zeros(growth, 1)];
    end
    
    % indices in dense vector
    indices = (nnz+1):(nnz+nnz_cs);
    % row indices in matrix (note we're building lower triangular, so shift these by i)
    r(indices) = i + nonzero_indices; 
    % this is correlation for ith variable
    c(indices) = i;
    
    s(indices) = nonzero_cs;
    
    nnz = nnz + nnz_cs;
end

if (nargout > 1)
    % knock off zero elements
    r(s == 0) = [];
    c(s == 0) = [];
    s(s == 0) = [];

else
    coef_mat = sparse(r(1:nnz), c(1:nnz), s(1:nnz), n_vars, n_vars, 2*nnz - n_vars);
    % isn't there just a property that a sparse matrix is symmetric? this seems
    % a bit wasteful
    s = coef_mat + coef_mat'; 
end
