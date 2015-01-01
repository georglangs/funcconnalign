function transform = align_points(S, T, varargin)
% aligns a source pointset to a target pointset
%
% [rotation, translation, S_mean, T_mean] = align_points(S, T, varargin)
%
% Arguments
% ---------
% S : n x p numeric
%   points to align from, each row is a point
% T : n x p numeric
%   points to align to, each row is a point
%
% Keyword Arguments
% -----------------
% weights : 0 <= n x 1 numeric
%   vector of weights, the ith weight corresponds to how much we care about aligning the ith source point (default = ones(n, 1))
% transform_type : {'rotation', 'orthonormal', 'rigid', 'scale'}
%   type of transformation, where:
%   'rotation' allows rotation only;
%   'orthonormal' allows rotation and flips;
%   'rigid' allows rotation and translation;
%   'scale' allows rotation, translation, and isotropic scale
%   (default = 'rigid')
%
% Returns
% -------
% transform : p x p double
%   rotation/orthonormal matrix from source to target points
% translation : 1 x p double
%   translation vector from source to target points

parser = inputParser();
parser.addRequired('S', @(x) validateattributes(x, {'double'}, {'2d', 'numeric'}));
parser.addRequired('T', @(x) validateattributes(x, {'double'}, {'2d', 'numeric'}));
parser.addParamValue('weights', [], @(x) validateattributes(x, {'double'}, {'vector', 'numeric', 'positive'}));
parser.addParamValue('transform_type', 'rigid', @(x) validatestring{'rotation', 'orthonormal', 'rigid', 'scale'}));
parser.parse(S, T, varargin{:});
inputs = parser.Results;

% further validation and initialization of inputs
[n, p] = size(T);

if n ~= size(S, 1)
    throw(fcalign.Exception('InvalidAttribute', 'Pointsets T and S must have the same number of points'));
end

if p ~= size(S, 2)
    throw(fcalign.Exception('InvalidAttribute', 'Pointsets T and S must have the same dimension'));
end

% what is this???
if (strcmp(transform_type, 'scale_orthonormal'))
    scale_correction = mean(sum(S.^2,2).^.5) / mean(sum(T.^2,2).^.5);
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
transform = zeros(p+1, p+1);
if strcmp(transform_type, 'orthonormal') || strcmp(transform_type, 'scale_orthonormal')
    transform(1:p, 1:p) = V*U';
else
    transform(1:p, 1:p) = V*diag([ones(1, p-1), det(V)*det(U)])*U';
end

% incorporate global scale correction
if strcmp(transform_type, 'scale_orthonormal')
    transform(1:p, 1:p) = scale_correction * transform(1:p, 1:p);
end

% calculate translation
if strcmp(transform_type, 'rigid')
    transform(1:p, p+1) = T_mean' - rotation*S_mean';
else
    transform(1:p, p+1) = zeros(p, 1);
end
