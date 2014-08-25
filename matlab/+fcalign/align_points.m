function [rotation, translation, S_mean, T_mean] = align_points(T, S, weights, transform_type)
% aligns a source pointset to a target pointset
%
% REQUIRED INPUTS
% T                     points to align from, each row is a point
%                       n x p double
% S                     points to align to, each row is a point
%                       n x p double
%
% OPTIONAL INPUTS
% weights               vector of weights, the ith weight corresponds to how much we care about aligning the ith source point
%                       n x 1 double (default = ones(n, 1))
% transform_type        type of transformation
%                       string form
%                           - 'rotation': determinant is 1
%                           - 'orthonormal': allow flips (determinant can be -1 or +1)
%                           - 'rigid': rotation and translation (default)
%                           - 'scale_orthonormal': incorporates global scale factor
%
% OUTPUTS
% rotation              rotation/orthonormal matrix from source to target points
%                       p x p double
% translation           translation vector from source to target points
%                       1 x p double
% S_mean                mean of source pointset
%                       1 x p double
% T_mean                mean of target pointset
%                       1 x p double
%
% Andy Sweet 11/02/2011

% parse inputs
parser = inputParser;
parser.addRequired('T', @(x) validateattributes(x, {'double'}, {'2d', 'real'}));
parser.addRequired('S', @(x) validateattributes(x, {'double'}, {'2d', 'real'}));
parser.addParameter('weights', [], @(x) validateattributes(x, {'double'}, {'vector', 'real', 'positive'}));
parser.addParameter('transform_type', 'rigid', @(x) validatestring{'rotation', 'orthonormal', 'rigid', 'scale_orthonormal'}));
parser.parse();
inputs = parser.Results;

% further validation and initialization of inputs
[n, p] = size(T);

if n ~= size(S, 1)
    error('fcalign:align_points:numberNotEqual', 'Pointsets T and S must have the same number of points');
end

if p ~= size(S, 2)
    error('fcalign:align_points:dimensionNotEqual', 'Pointsets T and S must have the same dimension');
end

% what is this???
if (strcmp(transform_type, 'scale_orthonormal'))
    scale_correction = mean(sum(S.^2,2).^.5) / mean(sum(T.^2,2).^.5)
    T = T * scale_correction;
end

% only center when there is translation
if strcmp(transform_type, 'rigid')
    S_mean = sum(bsxfun(@times, S, weights), 1) / sum(weights);
    T_mean = sum(bsxfun(@times, T, weights), 1) / sum(weights);

    S_centered = bsxfun(@minus, S, S_mean);
    T_centered = bsxfun(@minus, T, T_mean);
else
    S_mean = zeros(1, p);
    T_mean = zeros(1, p);

    S_centered = S;
    T_centered = T;
end

% calculate weight SVD
W = spdiags(weights, 0, n, n);
[U, Sigma, V] = svd(S_centered'*W*T_centered);

% form rotation matrix
if strcmp(transform_type, 'orthonormal') || strcmp(transform_type, 'scale_orthonormal')
    rotation = V*U';
else
    rotation = V*diag([ones(1, p-1), det(V)*det(U)])*U';
end

% incorporate global scale correction
if  (strcmp(transform_type, 'scale_orthonormal'))
    rotation = rotation * scale_correction;
end

% calculate translation
if strcmp(transform_type, 'rigid')
    translation = T_mean' - rotation*S_mean';
else
    translation = zeros(p, 1);
end
