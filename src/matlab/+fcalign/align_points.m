function transform = align_points(source, target, varargin)
    % Aligns a source pointset to a target pointset
    %
    % Arguments
    % ---------
    % source : s x n numeric
    %   points to align from, each row is a point
    % target : s x n numeric
    %   points to align to, each row is a point
    %
    % Keyword Arguments
    % -----------------
    % weights : 0 <= n x 1 numeric
    %   vector of weights, the sth weight corresponds to how much we care about aligning the sth source point (default = ones(s, 1))
    % transform_type : {'rotation', 'orthonormal', 'rigid', 'scale'}
    %   type of transformation (default = 'rigid') where:
    %   - 'rotation' allows rotation only;
    %   - 'orthonormal' allows rotation and flips;
    %   - 'rigid' allows rotation and translation;
    %   - 'scale' allows rotation, translation, and isotropic scale
    %
    % Returns
    % -------
    % transform : n x n double
    %   rotation/orthonormal matrix from source to target points
    % translation : 1 x n double
    %   translation vector from source to target points

    parser = inputParser();
    parser.addRequired('source', @(x) validateattributes(x, {'numeric'}, {'2d', 'real'}));
    parser.addRequired('target', @(x) validateattributes(x, {'numeric'}, {'2d', 'real'}));
    parser.addParamValue('weights', [], @(x) validateattributes(x, {'numeric'}, {'vector', 'real', 'positive'}));
    parser.addParamValue('transform_type', 'rigid', @(x) validatestring(x, {'rotation', 'orthonormal', 'rigid', 'scale'}));
    parser.parse(source, target, varargin{:});
    inputs = parser.Results;

    % further validation and initialization of inputs
    [s, n] = size(target);

    if s ~= size(source, 1)
        throw(fcalign.Exception('InvalidAttribute', 'Pointsets T and S must have the same number of points'));
    end

    if n ~= size(source, 2)
        throw(fcalign.Exception('InvalidAttribute', 'Pointsets T and S must have the same dimension'));
    end
    
    if isempty(inputs.weights)
        inputs.weights = ones(s, 1);
    end

    % what is this???
    if (strcmp(inputs.transform_type, 'scale_orthonormal'))
        scale_correction = mean(sum(source.^2,2).^.5) / mean(sum(target.^2,2).^.5);
        target = target * scale_correction;
    end

    % only center when there is translation
    if strcmp(inputs.transform_type, 'rigid')
        S_mean = sum(bsxfun(@times, source, inputs.weights), 1) / sum(inputs.weights);
        T_mean = sum(bsxfun(@times, target, inputs.weights), 1) / sum(inputs.weights);

        S_centered = bsxfun(@minus, source, S_mean);
        T_centered = bsxfun(@minus, target, T_mean);
    else
        S_mean = zeros(1, n);
        T_mean = zeros(1, n);

        S_centered = source;
        T_centered = target;
    end

    % calculate weight SVD
    W = spdiags(inputs.weights, 0, s, s);
    [U, Sigma, V] = svd(S_centered' * W * T_centered);

    % form rotation matrix
    transform = eye(n+1, n+1);
    if strcmp(inputs.transform_type, 'orthonormal') || ...
            strcmp(inputs.transform_type, 'scale_orthonormal')
        transform(1:n, 1:n) = V*U';
    else
        transform(1:n, 1:n) = V*diag([ones(1, n-1), det(V)*det(U)])*U';
    end

    % incorporate global scale correction
    if strcmp(inputs.transform_type, 'scale_orthonormal')
        transform(1:n, 1:n) = scale_correction * transform(1:n, 1:n);
    end

    % calculate translation
    if strcmp(inputs.transform_type, 'rigid')
        transform(1:n, n+1) = T_mean' - transform(1:n, 1:n) * S_mean';
    else
        transform(1:n, n+1) = zeros(n, 1);
    end

end