function [r ll] = calculate_posterior_prob_wrt_atlasclusters(x,clusters)
% function to calculate posterior probabilities for data points
% given clusters
%
% input:
% x .... points
% clusters ... cluster structure
% it is output of danials softKmeans.m
% 19.4.2011
% function derived from softKmeans.m
%%
listll = clusters.listll,%: [1x20 double]
m = clusters.m;%
p = clusters.p;%: [10x1 double]
sig = clusters.sig;%: 0.7706
r = clusters.r;%: [10x41501 double]
c= clusters.c;%: [1x41501 double]
ll= clusters.ll;%: [1x309 double]
k = size(m,1)
d = size(m,2)
if length(sig)==1
    sig = repmat(sig,k,d);
elseif min(size((sig))) ==1
    sig = repmat(sig,1,d);
end

x2 = x.^2;
m2 = m.^2;
bet = 1./sig;
n = size(x,1)
%clusters

% Compute negative of log likelihood matrix for each data point a cluster:
% nllmat : k x n

nllmat = 0.5*bet*x2'-(m.*bet)*x'+0.5*repmat( sum( m2.*bet - log(bet),2) ,1,n);
minNllmat = min(nllmat);

r = exp(-(nllmat-repmat(minNllmat,k,1)));



z = r'*p;
invz = 1./z; % (z+1e-15*~z);
invz(isinf(invz)) = 1e20;
r = r.*repmat(invz',k,1);

ll = -sum(log(invz))-sum(minNllmat);

%
if ~isempty(find(isnan(r(:))))
    keyboard
end
%}

end