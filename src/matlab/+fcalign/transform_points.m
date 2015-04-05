function Y = transform_points(X, T)
% Transforms points
%
% Arguments
% ---------
% X : n x s numeric
%   s n-dimensional points.
% T : (n+1) x (n+1)
%   Transform matrix.
%
% Returns
% -------
% Y : n x s numeric
%   Transformed points.

    parser = inputParser();
    parser.addRequired('X', @(x) validateattributes(x, {'numeric'}, {'2d', 'real'}));
    parser.addRequired('T', @(x) validateattributes(x, {'numeric'}, {'square', 'real'}));
    parser.parse(X, T);

    [n, s] = size(X);
    
    nT = size(T, 1);
    
    % check dimensions are consistent
    if (n + 1) ~= nT
        throw(fcalign.Exception('DimensionMismatch', [...
            'The dimension of the points is %u, but the dimension of the ', ...
            'transform is %u.'], n, nT ...
        ));
    end
            
    % form homogeneous coordinates
    X_h = [X; ones(1, s)];

    % transform homogeneous coordinates
    Y_h = T * X_h;

    % convert back to cartesian coordinates
    Y = bsxfun(@rdivide, Y_h(1:n, :), Y_h(n+1, :));
end
