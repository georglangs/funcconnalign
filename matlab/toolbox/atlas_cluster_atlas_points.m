function [atlas] = atlas_cluster_atlas_points(atlas,p)
% function to cluster atlas points
%
% georg langs
% 18.11.2011
%
%%
dp.cluster_method = 'soft_pcl'
dp.cluster_entire_atlas = true;
dp.cluster_individual_maps = false;
dp.n_clusters = 5;
dp.n_repeats = 20;
p = defaultParams(p,dp)

n_subjects = length(atlas.SubjectIndices);

%% Cluster entire Atlas

if (p.cluster_entire_atlas )
    if (strcmp(p.cluster_method, 'soft_j1'))
        clusters_struct = softKmeans(atlas.Gamma, p.n_clusters, p.n_repeats, 'j1', 'm');
    elseif (strcmp(p.cluster_method, 'soft_pcl'))
        clusters_struct = softKmeans(atlas.Gamma, p.n_clusters, p.n_repeats, 'pcl', 'm');
    elseif (strcmp(p.cluster_method, 'soft_pcldim'))
        clusters_struct = softKmeans(atlas.Gamma, p.n_clusters, p.n_repeats, 'pcldim', 'm');
    elseif (strcmp(p.cluster_method, 'soft_covars'))
        clusters_struct = softKmeans_fullvar(atlas.Gamma, p.n_clusters, [], p.n_repeats, 'j1', 'm');
    end
    
  %  eval(['atlas.', cluster_method, '.GroupClusterStructs{', num2str(n_clusters) '} = clusters_struct;']);
end

atlas.cluster_labels.atlas = clusters_struct.c;

if (p.cluster_individual_maps)
    %%
    for i=1:n_subjects
        
        subject_gamma = atlas.Gamma(atlas.SubjectIndices{i}, :);
        
        if (strcmp(p.cluster_method, 'soft_j1'))
            clusters_struct = softKmeans(subject_gamma, p.n_clusters, p.n_repeats, 'j1', 'm');
        elseif (strcmp(p.cluster_method, 'soft_pcl'))
            clusters_struct = softKmeans(subject_gamma, p.n_clusters, p.n_repeats, 'pcl', 'm');
        elseif (strcmp(p.cluster_method, 'soft_pcldim'))
            clusters_struct = softKmeans(subject_gamma, p.n_clusters, p.n_repeats, 'pcldim', 'm');
        elseif (strcmp(p.cluster_method, 'soft_covars'))
            clusters_struct = softKmeans_fullvar(subject_gamma, p.n_clusters, [], p.n_repeats, 'j1', 'm');
        end
        
    %    eval(['atlas.', cluster_method, '.SubjectClusterStructs{', num2str(i), ',' num2str(p.n_clusters) '} = clusters_struct;']);
        atlas.cluster_labels.individual(atlas.SubjectIndices{i},1) = clusters_struct.c;
    end
end


