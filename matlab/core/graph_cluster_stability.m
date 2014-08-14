function [rs, inst_rs, Rs] = graph_cluster_stability(W, Hs, ts, V, l)
% calculates stability of a set of clustered graphs over a range of transition times. Times must be increasing.

[n_cs, n_ts] = size(Hs);

n = size(W,1);

d = full(sum(W,2));
inv_sqrt_D = spdiags(d.^-0.5, 0, n, n);
sqrt_D = spdiags(d.^0.5, 0, n, n);

Vtsqrt_D = V'*sqrt_D;
inv_sqrt_DV = inv_sqrt_D*V;

pi = d' / sum(d);
Pi = spdiags(pi', 0, n, n);

% if there are fractional Markov times, check we have eigenvectors/values of adjoint
if ((ts(1) < 1) && (nargin < 5))
    error('Fractional Markov times require the eigenvectors/values of the adjoint of M.');
end

% compute stability for each number of clusters at each time step
Rs = cell(n_cs, n_ts);
rs = ones(n_cs, n_ts);
inst_rs = ones(n_cs, n_ts);
for i=1:n_cs
    
    tic;
    for j=1:n_ts

        H = Hs{i,j};
        c = size(H,2);

        % precompute some matrices and vectors we use a lot
        piH = pi*H; % 1 x c
        HtPi = H'*Pi; % c x n
    
        % compute t-step transition matrix product
        t = ts(j);      
        p = numel(l);
        L = spdiags(l.^t, 0, p, p);
        MH = inv_sqrt_DV*L*Vtsqrt_D*H;

        % clustered autocovariance matrix
        Rs{i,j} = HtPi*MH - piH'*piH;
        
        % then stability
        inst_rs(i,j) = sum(diag(Rs{i,j}));
        
        rs(i,j) = min([rs(i,1:(j-1)), inst_rs(i,j)]);
    end
    disp(['k=', num2str(c), ' took ', num2str(toc), ' s']);
end
