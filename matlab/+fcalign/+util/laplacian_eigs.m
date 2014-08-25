function [V, l] = laplacian_eigs(W, k, type)
% finds the eigendecomposition of the normalized Laplacian of a matrix
%
% [V, l] = diffusion_eigs(W, k)
%
% INPUTS
% W         affinity matrix
%           n x n double
% k         rank of eigendecomposition to find
%           1 <= integer <= n
% center    if true, centers laplacian before normalizing
%           logical
%
% OUTPUTS
% V         k principal eigenvectors
%           k x n double
% l         k principal eigenvalues
%           k x 1 double

n = size(W,1);

% compute laplacian and compute sums
if center
    W = W - diag(diag(W));
end
d = full(sum(W,2));
inv_sqrt_D = spdiags(d.^-0.5, 0, n, n);

% compute normalized laplacian
Ms = inv_sqrt_D*W*inv_sqrt_D;

% ensure symmetry (this is a hack to get round numerical precision lost above)
Ms = (Ms + Ms') / 2;

% get sorted eigen decomposition
[Vu, Lu] = eigs(Ms, k);
% to avoid forced symmetry use line below, but then MATLAB complains...
% [V,L] = sorted_eigs(Ms, k, 'lr');

lu = diag(Lu);
[l, l_inds] = sort(lu, descend);
V = Vu(:, l_inds);
