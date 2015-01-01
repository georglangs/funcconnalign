function Y = transform_points(X, T)
% transforms points
%
% Arguments
% ---------
% X : n x p numeric
%   n p-dimensional points
% T : (p+1) x (p+1)
%   transform matrix
%
% Returns
% -------
% Y : n x p numeric
%   transformed points
    parser = inputParser();
    parser.addRequired('X', @(x) validateattributes(x, {'numeric'}, {'matrix'}));
    parser.addRequired('T', @(x) validateattributes(x, {'numeric'}, {'square'}));
    parser.parse(X, T);

    [n, p] = size(X);

    % form homogeneous coordinates
    X_h = zeros(n, p+1);

    % transform homogeneous coordinates
    Y_h = T * X_h;

    % convert back to cartesian coordinates
    Y = bsxfun(@rdivide, Y_h(:, 1:p), Y_h(:, p+1));
end
