function gamma_new = diffusion_map_restriction(x_new, w_new, d, V, l)
% diffusion_map_restriction effectively interpolates in the diffusion map space to find a diffusion map coordinate for a new data vector
% 
% INPUTS
% x_new     new data vector
% w_new     affinities of new data vector with existing ones
% d         degrees of existing data
% V         eigenvectors of symmetric matrix adjoint to Markov chain matrix
% l         eigenvalues of symmetric matrix adjoint to Markov chain matrix
% 
% OUTPUTS
% gamma_new coordinates of new data vector in existing diffusion map
%
% See also diffusion_map_lifting
%

[n,k] = size(V);

[m,p] = size(x_new);

d_new = sum(w_new);

s_new = (d.^-0.5) .* (d_new.^-0.5) .* w_new;

v_new = zeros(k,1);

for i=1:k
    v_new(i) = sum(s_new.*V(:,i)) / l(i);
end

gamma_new = (d_new.^-0.5) .* v_new .* l;
