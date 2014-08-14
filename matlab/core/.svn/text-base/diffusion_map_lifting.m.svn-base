function x_new = diffusion_map_lifting(gamma_new, X, W, d, V, l, n_neighbors)
%  diffusion_map_lifting approximates the data vector that maps to the given diffusion map coordinate
% 
% INPUTS
% gamma_new diffusion map coordinate of new node
% X         existing data vectors
% W         affinities between existing data vectors
% d         degrees of existing nodes
% V         eigenvectors of symmetric matrix adjoint to Markov chain matrix
% l         eigenvalues of symmetric matrix adjoint to Markov chain matrix
% n_neighbors   number of neighbors to use for mean initialization
% 
% OUTPUTS
% x_new     approx data vector of new node
%
% See also diffusion_map_restriction
%

[n,k] = size(V);

inv_sqrt_D = spdiags(d.^-0.5, 0, n, n);
L = spdiags(l, 0, k, k);
Gamma = inv_sqrt_D*V*L;

% initialize solution with weighted mean of nearest neighbours
distances = distance(Gamma', gamma_new);
[sorted_distances, sorted_indices] = sort(distances, 'ascend');
nearest_distances = distances(1:n_neighbors);

init_x_new = sum(X(sorted_indices, :) ./ nearest_distances) / n_neighbors;

% then perform non-linear optimization around this solution
x_new = init_x_new;
