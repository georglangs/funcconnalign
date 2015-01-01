function map = create_diffusion_map(V, l, d, varargin)
% generates diffusion map coordinates from the eigendecomposition of the normalized laplacian
%
% Also normalizes the left and right eigenvectors so they have norm 1 w.r.t. the stationary distribution.
% Note that exactly 2 of p, t and delta must be provided (the other is automatically determined)
%
% Arguments
% ---------
% V : n x k double
%   eigenvectors of normalized laplacian
% l : 1 x k double
%   (1 - eigenvalues) of normalized laplacian
% d : n x 1 double
%   degrees of correlation matrix
%
% Keyword Arguments
% -----------------
% dim : 1 <= integer <= k
%   number of map dimensions (default = k)
% t : 1 <= double
%   diffusion time (default = 1)
% delta : double
%   threshold of ratio between extreme eigenvalues
%
% Returns
% -------
% map : struct
%   diffusion map struct containing:
%   Gamma - diffusion map coordinates (including redundant first coordinate);
%   determined_param - parameter left empty on input and determined from other 2;
%   Psi - right eigenvectors of Markov matrix;
%   Phi - left eigenvectors of Markov matrix;
parser = inputParser();
parser.addRequired('V', @(x) validateattributes(x, {'double'}, {'2d'}));
parser.addRequired('l', @(x) validateattributes(x, {'double'}, {'vector'}));
parser.addRequired('d', @(x) validateattributes(x, {'double'}, {'vector'}));
parser.addParamValue('dim', [], @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
parser.addParamValue('t', [], @(x) validateattributes(x, {'double'}, {'scalar', 'positive'}));
parser.addParamValue('delta', [], @(x) validateattributes(x, {'double'}, {'scalar'}));
parser.parse(V, l, d, varargin{:});
inputs = parser.Results;

% further process inputs
[n, k] = size(V);

dim = inputs.dim;
t = inputs.t;
delta = inputs.delta;

if length(l) ~= k
    throw(fcalign.Exception('InvalidAttribute', 'length of l must be same as number of columns in V'));
end

if length(d) ~= n
    throw(fcalign.Exception('InvalidAttribute', 'length of d must be same as number of rows in V'));
end

if (isempty(dim) + isempty(t) + isempty(delta)) ~= 1
    throw(fcalign.Exception('TooManyArguments', 'exactly 2 of dim, t and delta must be provided'));
end

% determine one of the parameters
if isempty(dim)
    % use diffusion time and delta to determine dimensionality
    for dim = 1:numel(l)
        delta_try = (l(dim) / l(2))^t;
        if ( delta_try <= delta )
            break;
        end
    end

    map.diffusion_time = t;
    map.delta = delta;
    map.map_dimension = dim;
    map.determined_param = 'dim';

elseif isempty(t)
    % use dimensionality and delta to determine diffusion time
    t = log(delta) / (log(l(p.dim) / l(2)));

    map.diffusion_time = t;
    map.delta = delta;
    map.map_dimension = dim;
    map.determined_param = 'diffusion_time';

elseif isempty(delta)
    % use diffusion time and dimensionality to determine delta
    delta = (l(dim) / l(2))^t;

    map.diffusion_time = t;
    map.delta = delta;
    map.map_dimension = dim;
    map.determined_param = 'delta';
end

% this ensures that norm wrt stationary distribution (phi_0) is 1
norm_factor = sqrt(sum(d));

% construct left and right eigenvectors of stochastic matrix from
inv_sqrt_D = spdiags(d.^-0.5, 0, n, n);
sqrt_D = spdiags(d.^0.5, 0, n, n);
map.Phi = sqrt_D*V(:, 1:dim) / norm_factor;
map.Psi = norm_factor * inv_sqrt_D*V(:, 1:dim);

% first right eigenvector should be constant positive, so if it isn't positive, we need to flip all eigenvectors
if map.Psi(1, 1) < 0
    map.Phi = map.Phi * -1;
    map.Psi = map.Psi * -1;
end

% construct diffusion map coordinates
Lt = spdiags(l(1:dim).^t, 0, dim, dim);
map.Gamma = map.Psi*Lt;
