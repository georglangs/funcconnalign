function funcmap_visualize2(map_points, subplot_dims, colors, vis_sample)
% visualizes a functional map over some permutations of its dimensions
% points is n x d x k:, n is number of points in cluster, d is dimension of points, and k is cluster index

if (nargin < 4)
    vis_sample = 1;
end

m = subplot_dims(1);
n = subplot_dims(2);
n_clusters = numel(map_points);

for i=1:m
    
    i_max = 0;
    for k=1:n_clusters
        i_max_cluster = max(abs(map_points{k}(:,i)));
        if (i_max_cluster > i_max)
            i_max = i_max_cluster;
        end
    end
    
    for j=1:n
        
        index = (i-1)*n + j;
        
        j_max = 0;
        for k=1:n_clusters
            j_max_cluster = max(abs(map_points{k}(:,i+j)));
            if (j_max_cluster > j_max)
                j_max = j_max_cluster;
            end
        end
        
        subplot(m,n,index), hold on;
        for k=1:n_clusters
            
            cluster_points = map_points{k};
            
            plot(cluster_points(1:vis_sample:end, i), cluster_points(1:vis_sample:end, i+j), '.', 'color', colors(k,:));
            
        end
        
        title([num2str(i), ', ', num2str(i+j)]); 
        xlim([-i_max,i_max]), ylim([-j_max,j_max]);
        axis equal;
    end
end
