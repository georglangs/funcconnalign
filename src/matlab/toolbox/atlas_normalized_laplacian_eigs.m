function [V, l] = atlas_normalized_laplacian_eigs(W, p)
% ATLAS_NORMALIZED_LAPLACIAN_EIGS find the p-rank eigendecomposition a matrix that is essentially equivalent to the normalized Laplacian
%
% INPUTS
% W     n x n: affinity matrix
% p     rank of eigendecomposition to find
% 
% OUTPUTS
% V     p x 1: prinicpal eigenvectors of normalized Laplacian
% l     p x 1: 1 - principal eigenvalues of normalized Laplacian

n = size(W,1);


d = full(sum(W,2));
inv_sqrt_D = spdiags(d.^-0.5, 0, n, n);

% normalized Laplacian is inv_sqrt_D*(D-W)*inv_sqrt_D, but here we ignore the identity matrix 

Ms = inv_sqrt_D*W*inv_sqrt_D;

% ensure symmetry (this is a hack to get round numerical precision lost above)
Ms = (Ms + Ms') / 2;
[V,L] = sorted_eigs(Ms, [], p);

% to avoid use line below, but then MATLAB complains...
% [V,L] = sorted_eigs(Ms, [], p, 'lr');

l = diag(L);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V,L] = sorted_eigs(A, B, k, sigma, eig_opts, order)
% SORTED_EIGS computes the k-rank generalized eigendecomposition of A wrt B: i.e. A*V = B*V*L
%
% INPUTS
% A         n x n: sparse square matrix 
% B         n x n: SPD matrix
% k         rank of eigendecomposition
% sigma     see eigs: sigma
% eig_opts  see eigs: opts
% order     sort order
% 
% OUTPUTS
% V         n x k: generalized eigenvectors of A wrt B
% L         k x k: generalized eigenvalues of A wrt B
%
% written by Andy Sweet, 16/05/2011

gen = true;
if ~exist('B', 'var') || isempty(B)
    gen = false;
end

if ~exist('order', 'var') || isempty(order)
    order = 'descend';
end

if ~exist('sigma', 'var') || isempty(sigma)
    sigma = 'la';
end

if ~exist('eig_opts', 'var') || isempty(eig_opts)
    eig_opts.issym = true;
    eig_opts.isreal = true;
end

% unsorted eigendecomp
if gen
    [Vu, Lu] = eigs(A, B, k, sigma, eig_opts);
else
    [Vu, Lu] = eigs(A, k, sigma, eig_opts);
end

lu = diag(Lu);

% sort in specified order
[l, l_inds] = sort(lu, order);
L = diag(l);
V = Vu(:, l_inds);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V,L] = sorted_eig(A, B, k, sigma, eig_opts, order)
% SORTED_EIGS computes the k-rank generalized eigendecomposition of A wrt B: i.e. A*V = B*V*L
%
% INPUTS
% A         n x n: sparse square matrix 
% B         n x n: SPD matrix
% k         rank of eigendecomposition
% sigma     see eigs: sigma
% eig_opts  see eigs: opts
% order     sort order
% 
% OUTPUTS
% V         n x k: generalized eigenvectors of A wrt B
% L         k x k: generalized eigenvalues of A wrt B
%
% written by Andy Sweet, 16/05/2011

gen = true;
if ~exist('B', 'var') || isempty(B)
    gen = false;
end

if ~exist('order', 'var') || isempty(order)
    order = 'descend';
end

if ~exist('sigma', 'var') || isempty(sigma)
    sigma = 'la';
end

if ~exist('eig_opts', 'var') || isempty(eig_opts)
    eig_opts.issym = true;
    eig_opts.isreal = true;
end

% unsorted eigendecomp
if gen
    [Vu, Lu] = eig(full(A));
else
    [Vu, Lu] = eig(full(A));
end
Vu = Vu(:,end-k+1:end);
Lu = Lu(end-k+1:end,end-k+1:end);
lu = diag(Lu);

lu = lu(end:-1:1);
Vu = Vu(:,end:-1:1);

% sort in specified order
[l, l_inds] = sort(lu, order);
L = diag(l);
V = Vu(:, l_inds);
