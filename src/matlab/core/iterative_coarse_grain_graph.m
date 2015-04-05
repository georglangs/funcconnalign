function [merged_labels, merged_P, merged_W] = iterative_coarse_grain_graph(P, phi0, labels)
% iterative_coarse_grain_graph iteratively merge coarse grain graph until each cluster prefers to self transition
%
% INPUTS
% P         Markov chain matrix of original graph (n x n)
% phi0      stationary distribution of original graph (n x 1)
% labels    initial subgraph labels for each node  (n x 1)
%
% OUTPUTS
% merged_labels     merged subgraph labels
% merged_P          Markov chain matrix of merged coarse grain graph
% merged_W          weight/similarity matrix of merged coarse grain graph

merged_labels = labels;
converged = 0;
while ~converged
    
    % get new coarse grain graph based on current cluster labels
    [merged_P, merged_W] = coarse_grain_graph(P, phi0, merged_labels);
    
    % get max transition probs for each cluster
    [max_ps, max_cols] = max(merged_P, [], 2);
    
    self_ps = diag(merged_P);
    
    % find those clusters who prefer to transition to another (rather than self transition)
    n_clusters_merged = max(merged_labels);
    non_self_transitioners = ((1:n_clusters_merged)' ~= max_cols);
    
    % and change the one that has the strongest aversion to itself
    if (sum(non_self_transitioners) > 0)
        
        ratio_ps = self_ps./max_ps;
        [min_ratio, cluster_to_merge] = min(ratio_ps);
        parent_cluster = max_cols(cluster_to_merge);
        
        merged_labels(merged_labels == cluster_to_merge) = parent_cluster;
        
        % assign largest remaining cluster index to merged cluster
        if (n_clusters_merged ~= cluster_to_merge)
            merged_labels(merged_labels == n_clusters_merged) = cluster_to_merge;
        end
        
    else
    % or finish
        converged = 1;
    end
end
