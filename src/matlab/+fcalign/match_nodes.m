function [matched_indices, matched_indices_value, visu_vertices, visu_indices] = match_nodes(X, Y, thresh, thresh2,max_nr_pairs)
% ATLAS_SPARSE_CROSS_DISTANCE creates a sparse square euclidean distance matrix given two data matrices X and Y and an upper bound thresh
% 
% INPUTS
% X                 n x p: data matrix where each row is an observation (required)
% Y                 n x p: data matrix where each row is an observation (required)
% thresh            upper threshold on distances or significance values (required)
% thresh_type       what to threshold
%                       'distance': Euclidean distance (default)
%                       'significance': Bonferroni corrected p-value
% nzmaxi             the maximum number of non-zero values we expect to get
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
%%
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

if ~exist('thresh_type', 'var') || isempty(thresh_type)
    thresh_type = 'distance';
end

if ~exist('growth_rate', 'var') || isempty(growth_rate)
    growth_rate = 0.25;
end

% if we're not given the expected max number of nonzeros, then assume sparsity won't exceed 1/50
if ~exist('nzmaxi', 'var') || isempty(nzmaxi)
    nzmaxi = ceil((p^2)/50);
    if verbose
        disp(['using auto nzmaxi ', num2str(nzmaxi)]);
    end
end

r = zeros(nzmaxi, 1);
c = zeros(nzmaxi, 1);
s = zeros(nzmaxi, 1);

% we grow these these vectors if they're full at a rate proportional to the original max
growth = int32(growth_rate*nzmaxi);

% find distance by variable
nnz = 0;
for i=1:p_x
    
    if verbose
        if (~mod(i,500))
            disp(['calculated for ', num2str(i), ' variables']);
        end
    end
    
    % find distances between this x variable and all ys
    ds = sqrt(sum((repmat(X(:,i), [1, p_y]) - Y).^2));
    
    % find the values that should be nonzero and their inds
    if (strcmp(thresh_type, 'significance'))

        % TODO: use variance estimate to find significant distances?
        % n_dofs = n - 2;
        % n_comparisons = p_x*p_y;
        % 
        % ts = cs .* (n_dofs ./ (1 - cs.^2)).^0.5;
        % 
        % t_thresh = tinv(1 - thresh/n_comparisons, n_dofs);
        % nonzero_inds = find(ts > t_thresh);
        error('Significance test not implement yet ...');
        
    else
        nonzero_inds = find(ds < thresh);
    end
    
    nonzero_ds = ds(nonzero_inds);
    
    % if adding these takes us over max, then add some more space
    % (proportion to original max)
    nnz_ds = length(nonzero_ds);
    if ( nnz + nnz_ds > length(s) )
        r = [r; zeros(growth, 1)];
        c = [c; zeros(growth, 1)];
        s = [s; zeros(growth, 1)];
    end
    
    % inds in dense vector
    inds = (nnz+1):(nnz+nnz_ds);
    % row inds in matrix (note we're building lower triangular, so shift these by i)
    r(inds) = i; 
    % this is correlation for ith variable
    c(inds) = nonzero_inds;
    
    s(inds) = nonzero_ds;
    
    nnz = nnz + nnz_ds;
end

%if (nargout > 1)
    % knock off zero elements
    r(s == 0) = [];
    c(s == 0) = [];
    s(s == 0) = [];
%else
%    s = sparse(r(1:nnz), c(1:nnz), s(1:nnz), p_x, p_y, nnz);
%end
disp('...')
%% Perform matching:
disp('Starting to match ....')
   got_the_node = zeros(size(X,2),2);
counter = 1;
counter_matched = 1;
matches_temp = [r c s]; % row index, column index, value
matches_temp = sortrows(matches_temp,3);
matched_indices = [];
if not(exist('max_nr_pairs','var'))
    max_nr_pairs = size(X,2)+1;
end
%h = waitbar(0)
while min(got_the_node) == 0 & counter_matched < size(X,2) & counter < size(matches_temp,1) & matches_temp(counter,3) < thresh2 & counter_matched < max_nr_pairs
    % spy(got_the_node)
    % drawnow
    
    if counter > 1
        if (not(any(matched_indices(:,1) == matches_temp(counter,1))) & not(any(matched_indices(:,2) == matches_temp(counter,2))))
            % if the indices have not been used up
            matched_indices(counter_matched,1:2) = matches_temp(counter,1:2);
            matched_indices_value(counter_matched,1) = matches_temp(counter,3);
            got_the_node(matched_indices(counter_matched,1),1) = 1;
            got_the_node(matched_indices(counter_matched,2),2) = 1;
            %disp(num2str(matched_indices(counter_matched,1:2)));
            counter_matched
            counter_matched = counter_matched + 1;

        end
    else
        matched_indices(counter_matched,1:2) = matches_temp(counter,1:2);
        matched_indices_value(counter_matched,1) = matches_temp(counter,3);
        got_the_node(matched_indices(counter_matched,1),1) = 1;
        got_the_node(matched_indices(counter_matched,2),2) = 1;
        %disp(num2str(matched_indices(counter_matched,1:2)));
        counter_matched
        counter_matched = counter_matched + 1;
    end
    counter = counter + 1;
end
if (counter_matched == size(X,2))
    disp('.... nodes matched successfully')
else
    disp('.... finished but not all nodes matched')
    
end

%{
visu_vertices = [];
visu_indices = [1;2];
n_nodes = size(X,1);
for i = 1:size(matched_indices,1)
    visu_vertices = [visu_vertices; X(:,matched_indices(i,1))'; Y(:,matched_indices(i,2))'];
   visu_indices(i,1) = [2*i-1];% [visu_indices ; matched_indices(i,2); n_nodes+matched_indices(i,1) ];
end
%}



