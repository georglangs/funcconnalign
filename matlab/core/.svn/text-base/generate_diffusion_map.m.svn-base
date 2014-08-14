function [Gamma, determined_param, Psi, Phi] = generate_diffusion_map(V, l, d, p, t, delta)
% generates diffusion map coordinates from the eigendecomposition of the normalized laplacian. Also normalizes the left and right eigenvectors so they have norm 1 w.r.t. the stationary distribution.
%
% INPUTS
% V     eigenvectors of normalized laplacian
% l     eigenvalues of normalized laplacian
% d     degrees of affinity matrix
% p     number of dimensions to take
% t     diffusion time
% delta threshold of ratio between extreme eigenvalues
%
% OUTPUTS
% Gamma diffusion map coordinates (including redundant first coordinate)
% determined_param  parameter left empty on input and determined from other 2
% Psi   right eigenvectors of Markov matrix
% Phi   left eigenvectors of Markov matrix
%

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
    disp('Either dimensionality, diffusion time or delta threshold must be left unspecified.')
end

n = numel(d);
inv_sqrt_D = spdiags(d.^-0.5, 0, n, n);
sqrt_D = spdiags(d.^0.5, 0, n, n);
norm_factor = sqrt(sum(d));

Phi = sqrt_D*V(:, 1:p) / norm_factor;
Psi = norm_factor * inv_sqrt_D*V(:, 1:p);

% negate normalization factor if first eigenvector is negative
if (Psi(1,1) < 0)
    Phi = Phi * -1;
    Psi = Psi * -1;
end

Lt = spdiags(l(1:p).^t, 0, p, p);
Gamma = Psi*Lt;
