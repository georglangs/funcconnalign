function [V, l] = normalized_laplacian_eigs(W, p)
% generate_diffusion_map creates a p-dimensional diffusion map from an affinity matrix
%
% INPUTS
% W     affinity matrix
% p     dimension of diffusion map (how many eigenvectors to find)
% 
% OUTPUTS
% V     p prinicpal eigenvectors
% l     p principal eigenvalues

n = size(W,1);

d = full(sum(W,2));
inv_sqrt_D = spdiags(d.^-0.5, 0, n, n);

Ms = inv_sqrt_D*W*inv_sqrt_D;

eig_opts.issym = 1;
[V,L] = eigs(Ms, p, 'lr', eig_opts);
[l, sorted_indices] = sort(diag(L),'descend');

V = V(:, sorted_indices);
