function [rotation, translation, S_mean, T_mean] = atlas_weighted_align_shapes_nd(T, S, weights, transform_type)
% ATLAS_WEIGHTED_ALIGN_SHAPES_ND aligns a source pointset to a target pointset 
%
% INPUTS
% T                     n x p: points to align from, each row is a point
% S                     n x p: points to align to, each row is a point
% weights               n x 1: vector of weights, the ith weight corresponds to how much we care about aligning the ith source point
% transform_type        type of transformation
%                           - 'rotation': determinant is 1
%                           - 'orthonormal': allow flips (determinant can be -1 or +1)
%                           - 'rigid': rotation and translation
%
% OUTPUTS
% rotation              p x p: rotation/orthonormal matrix from source to target points
% translation           p x 1: translation vector from source to target points
% S_mean                p x 1: mean of source pointset
% T_mean                p x 1: mean of target pointset
%
% Andy Sweet 11/02/2011

[n_t, p_t] = size(T);
[n_s, p_s] = size(S);

if (p_t ~= p_s)
    error('dimensions of pointsets not equal')
end
p = p_s;

if n_s ~= n_t
    error('cardinality of pointsets not equal')
end
n = n_s;

if ~exist('weights', 'var') || isempty(weights)
    weights = ones(n,1);
end

if  (strcmp(transform_type, 'scale_orthonormal'))
 scale_correction = mean(sum(S.^2,2).^.5) / mean(sum(T.^2,2).^.5)
 T = T * scale_correction;
end

if (strcmp(transform_type, 'rotation') || strcmp(transform_type, 'orthonormal') || strcmp(transform_type, 'scale_orthonormal'))

    S_mean = zeros(1, p);
    T_mean = zeros(1, p);
    
    S_centered = S;
    T_centered = T;    
else
    S_mean = sum(bsxfun(@times, S, weights), 1) / sum(weights);
    T_mean = sum(bsxfun(@times, T, weights), 1) / sum(weights);
    
    S_centered = bsxfun(@minus, S, S_mean);    
    T_centered = bsxfun(@minus, T, T_mean);
end

W = spdiags(weights, 0, n, n);
[U,Sigma,V] = svd(S_centered'*W*T_centered);

if (strcmp(transform_type, 'orthonormal') || strcmp(transform_type, 'scale_orthonormal'))
    rotation = V*U';
else
    rotation = V*diag([ones(1, p-1), det(V)*det(U)])*U';
end

if  (strcmp(transform_type, 'scale_orthonormal'))
 rotation = rotation * scale_correction;
end

if (strcmp(transform_type, 'rotation') || strcmp(transform_type, 'orthonormal') || strcmp(transform_type, 'scale_orthonormal'))
    translation = zeros(p,1);
else
    translation = T_mean' - rotation*S_mean';
end
