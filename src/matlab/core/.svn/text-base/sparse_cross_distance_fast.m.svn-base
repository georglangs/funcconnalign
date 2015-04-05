function [s,r,c] = sparse_cross_distance_fast(X, Y, thresh, thresh_type, nzmax)
% sparse_corrcoef creates a sparse square euclidean distance matrix given a data matrix X and a distance upper bound
% 
% INPUTS
% X         data matrix where each row is an observation
% Y         data matrix where each row is an observation
% thresh    lower bound of correlation coefficients to keep (to achieve sparsity)
% nzmax     the maximum number of non-zero values we expect to get (optional)
%
% OUTPUTS
% s         sparse distance matrix (or, if 3 outputs specfied, the distances in a dense vector)
% r         row indices of sparse entries
% c         column indices of sparse entries
%

[n_obs_x, n_vars_x] = size(X);
[n_obs_y, n_vars_y] = size(Y);

if ( n_obs_x ~= n_obs_y )
    error(['Number of observations in X, ', num2str(n_obs_x), ', does not equal number of observations in Y, ', num2str(n_obs_y), '.']);
end
n_obs = n_obs_x;
n_vars = (n_vars_x + n_vars_y)/2;

if (nargin > 5)
    disp(['using nzmax ', num2str(nzmax)]);
else
    % if we're not given the expected max number of nonzeros, then assume
    % this has ratio 1/50
    nzmax = ceil((n_vars^2)/50);
    disp(['using auto nzmax ', num2str(nzmax)]);
end

if nargin < 4
    thresh_type = 'distance';
end

r = zeros(nzmax, 1);
c = zeros(nzmax, 1);
s = zeros(nzmax, 1);

% we grow these these vectors if they're full at a rate proportional to the original max
growth_rate = 0.25;
growth = int32(growth_rate*nzmax);

% find correlations by variable
nnz = 0;
for i=1:n_vars_x
    
    % find distances between this variable and all
    ds = sqrt(sum((repmat(X(:,i), [1, n_vars_y]) - Y).^2));
    
    % find the values that should be nonzero and their indices
    if (strcmp(thresh_type, 'significance'))

        % TODO: use variance estimate to find significant distances?
        % n_dofs = n_obs - 2;
        % n_comparisons = n_vars_x*n_vars_y;
        % 
        % ts = cs .* (n_dofs ./ (1 - cs.^2)).^0.5;
        % 
        % t_thresh = tinv(1 - thresh/n_comparisons, n_dofs);
        % nonzero_indices = find(ts > t_thresh);
        
    else
        nonzero_indices = find(ds < thresh);
    end
    
    nonzero_ds = ds(nonzero_indices);
    
    % if adding these takes us over max, then add some more space
    % (proportion to original max)
    nnz_ds = length(nonzero_ds);
    if ( nnz + nnz_ds > length(s) )
        r = [r; zeros(growth, 1)];
        c = [c; zeros(growth, 1)];
        s = [s; zeros(growth, 1)];
    end
    
    % indices in dense vector
    indices = (nnz+1):(nnz+nnz_ds);
    % row indices in matrix (note we're building lower triangular, so shift these by i)
    r(indices) = i; 
    % this is correlation for ith variable
    c(indices) = nonzero_indices;
    
    s(indices) = nonzero_ds;
    
    nnz = nnz + nnz_ds;
end

if (nargout > 1)
    % knock off zero elements
    r(s == 0) = [];
    c(s == 0) = [];
    s(s == 0) = [];
else
    s = sparse(r(1:nnz), c(1:nnz), s(1:nnz), n_vars_x, n_vars_y, nnz);
end
