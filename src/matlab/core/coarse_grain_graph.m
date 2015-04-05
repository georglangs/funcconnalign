function [Pk,Wk,phi0k] = coarse_grain_graph(P, phi0, labels)
% coarse_grain_graph creates the markov chain matrix Pk and weight matrix Wk, of a graph with markov chain matrix P and stationary distribution phi0 partitioned into k subgraphs (clusters)
%
% INPUTS
% P         Markov chain matrix of original graph (n x n)
% phi0      stationary distribution of original graph (n x 1)
% labels    subgraph labels for each node  (n x 1)
%
% OUTPUTS
% Pk        Markov chain matrix of coarse grain graph
% Wk        weight/similarity matrix of coarse grain graph
%

% find number of subgraphs/clusters
k = max(labels);

% calculate weights between each subgraph
Wk = zeros(k,k);

phi0k = zeros(k,1);

for i=1:k
    cluster_is = (labels == i);
    phi0k(i) = sum(phi0(cluster_is));
    for j=1:k
        cluster_js = (labels == j);
        Wk(i,j) = phi0(cluster_is)'*full(sum(P(cluster_is, cluster_js),2));
    end
end

% generate markov chain matrix from weight matrix
Pk = bsxfun(@rdivide, Wk, sum(Wk,2));
