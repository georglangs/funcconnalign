function [rotation, translation, source_mean, target_mean] = weightedAlignShapesND(target_points, source_points, weights, transform_type)
% rigidly aligns a source pointset to a target pointset 
%
% INPUTS
% target_points         points to align to, each row is a point index, each column a dimension
% source_points         points to align, each row is a point index, each column a dimension
% point_pair_weights    vector of weights, the ith weight corresponds to how much we care about aligning the ith source point?
% transform_type        type of transformation
%                           - 'rotation'
%                           - 'orthonormal'
%                           - 'rigid'
%
% OUTPUTS
% aligned_source_points
% rotation              rotation matrix from source to target points
% source_mean           mean of source pointset
% target_mean           mean of target pointset
%
% Andy Sweet 11/02/2011

[n_target_points, n_target_dims] = size(target_points);
[n_source_points, n_source_dims] = size(source_points);

if n_source_dims ~= n_target_dims
    disp('dimensions of pointsets not equal')
    return;
end
n_dims = n_source_dims;

if n_source_points ~= n_target_points
    disp('cardinality of pointsets not equal')
    return;
end
n_points = n_source_points;

if isempty(weights)
    weights = ones(n_points, 1);
end

if (strcmp(transform_type, 'rotation') || strcmp(transform_type, 'orthonormal'))

    source_mean = zeros(1, n_dims);
    target_mean = zeros(1, n_dims);
    
    centered_target_points = target_points;
    centered_source_points = source_points;    
else
    source_mean = sum(bsxfun(@times, source_points, weights), 1) / sum(weights);
    target_mean = sum(bsxfun(@times, target_points, weights), 1) / sum(weights);
    
    centered_target_points = bsxfun(@minus, target_points, source_mean);
    centered_source_points = bsxfun(@minus, source_points, target_mean);
end

weight_matrix = spdiags(weights, 0, n_points, n_points);
[U,S,V] = svd(centered_source_points'*weight_matrix*centered_target_points);

if (strcmp(transform_type, 'orthonormal'))
    rotation = V*U';
else
    rotation = V*diag([ones(1, n_dims-1), det(V)*det(U)])*U';
end

if (strcmp(transform_type, 'rotation') || strcmp(transform_type, 'orthonormal'))
    translation = zeros(n_dims,1);
else
    translation = target_mean' - rotation*source_mean';
end
