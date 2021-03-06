function [Gamma, determined_param, Psi, Phi] = atlas_generate_diffusion_map(V, l, d, p, t, delta)
% ATLAS_GENERATE_DIFFUSION_MAP generates diffusion map coordinates from the eigendecomposition of the normalized laplacian. Also normalizes the left and right eigenvectors so they have norm 1 w.r.t. the stationary distribution.
%
% INPUTS
% V         eigenvectors of normalized laplacian
% l         (1 - eigenvalues) of normalized laplacian
% d         degrees of affinity matrix
% p         number of dimensions
% t         diffusion time
% delta     threshold of ratio between extreme eigenvalues
%
% OUTPUTS
% Gamma     diffusion map coordinates (including redundant first coordinate)
% determined_param  parameter left empty on input and determined from other 2
% Psi       right eigenvectors of Markov matrix
% Phi       left eigenvectors of Markov matrix
%
% Andy Sweet 12/02/2011	 

% check which parameters determine
if (isempty(p))
    % use diffusion time and delta to determine dimensionality
    for p=1:numel(l)
        delta_try = (l(p) / l(2))^t;
        if ( delta_try <= delta )
            break;
        end
    end
    
    determined_param = p;
    
elseif (isempty(t))
    % use dimensionality and delta to determine diffusion time
    t = log(delta) / (log(l(p) / l(2)));
    
    determined_param = t;
    
elseif (isempty(delta))
    % use diffusion time and dimensionality to determine delta
    delta = (l(p) / l(2))^t;
    
    determined_param = delta;
else
    error('Either dimensionality, diffusion time or delta threshold must be left unspecified.')
end

n = numel(d);

% this ensures that norm wrt stationary distribution (phi_0) is 1
norm_factor = sqrt(sum(d));

% construct left and right eigenvectors of stochastic matrix from 
inv_sqrt_D = spdiags(d.^-0.5, 0, n, n);
sqrt_D = spdiags(d.^0.5, 0, n, n);
Phi = sqrt_D*V(:, 1:p) / norm_factor;
Psi = norm_factor * inv_sqrt_D*V(:, 1:p);

% first right eigenvector should be constant positive, so if it isn't positive, we need to flip all eigenvectors
if (Psi(1,1) < 0)
    Phi = Phi * -1;
    Psi = Psi * -1;
end

% construct diffusion map coordinates
Lt = spdiags(l(1:p).^t, 0, p, p);
Gamma = Psi*Lt;
