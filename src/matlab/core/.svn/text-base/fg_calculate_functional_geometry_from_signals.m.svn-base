function [map] = fg_calculate_functional_geometry_from_signals(data,p)
%
% function to calculate the functional geometry from a set of signals
% 
% INPUT: 
% S.BOLD ... signal matrix, one row= one signal, one column = one timepoint
% S.coords ... (optional) coordinates of the signal sources
% p ... parameter structure
%
% Functional Geometry Toolbox V0.1
% 14.2.2011
% langs@csail.mit.edu
% sweet@csail.mit.edu
%

if nargin == 0
    % generation of demo data if no inputs are provided
    T = linspace(0,2*pi,100)';
    S1 = 2*sin(T);
    S2 = (S1 + 1*randn(100,1));
    S3 = S2 + .6*randn(100,1);
    Noise = .6;
    RegionSize = 120;
    BGSize = 40;
    data.BOLD = [ ...
        .1*randn(BGSize,100); ...
        repmat(S1',RegionSize,1); ...
        repmat(S2',RegionSize,1); ...
        repmat(S3',RegionSize,1);
        ];
    data.BOLD = data.BOLD + Noise*randn(size(data.BOLD));
    plot(data.BOLD')
end

%% Setting default parameters:
dp.do_save_memory = 0;
dp.sampling_type = 'corr_sum';
dp.thresh_type =  'correlation';
dp.correlation_thresh = .6;
dp.significance_thresh = .05;
% diffusion map parameters
dp.n_dims = 150;
dp.epsilon = 0.65;
dp.correlation_thresh = 0.6;
dp.min_corr_sum = 100;
dp.min_density = -1; % if negative then automatically choose better than chance
dp.significance_thresh = 0.01;
dp.keep_negative_correlations = 0;
dp.alpha = 0;
dp.do_save_memory = 0;
dp.thresh_type = 'correlation';
% thresh_type = 'significance';

p.sampling_type = 'corr_sum'; 
% sampling_type = 'density';

p = defaultParams(p,dp);
%% generating sparse correlation coefficient matrix
data.NonZeroIndices = find(mean(data.BOLD,2)~=0);

disp('Generating correlation matrix ...');
tic;
if p.do_save_memory
    
    if (strcmp(p.thresh_type, 'significance'))
        GCT = sparse_corrcoef_fast(data.BOLD(data.NonZeroIndices,:)', p.significance_thresh, p.thresh_type, p.keep_negative_correlations);
        
    elseif (strcmp(p.thresh_type, 'correlation'))
        GCT = sparse_corrcoef_fast(data.BOLD(data.NonZeroIndices,:)', p.correlation_thresh, p.thresh_type, p.keep_negative_correlations);
    else
        disp('No other threshold types available yet. ')
        return;
    end
    
else
    GC = corrcoef(data.BOLD(data.NonZeroIndices,:)');
    
    if (strcmp(p.thresh_type, 'correlation'))
        GCT = sparse(double(GC.*(GC > p.correlation_thresh)));
    else
        disp('No other threshold types available yet. ')
        return;
    end
    
    GCT(isnan(GCT)) = 0;
    clear GC
end
disp(['... done in ', num2str(toc), 's']);

%% calculating the degree of the nodes, and discarding nodes with degree below threshold
disp('Generating and sampling affinity matrix ...'); tic;
n_nodes = numel(data.NonZeroIndices);

W = sparse(n_nodes, n_nodes);
W(GCT~=0) = exp(GCT(GCT~=0) / p.epsilon);

if (strcmp(p.sampling_type, 'corr_sum'))
    
    corr_sums = full(sum(GCT));
    nodes_kept = find(corr_sums > p.min_corr_sum);
    
elseif (strcmp(p.sampling_type, 'density'))
    
    d_orig = full(sum(W,2));
    pi_orig = d_orig / sum(d_orig);
    
    if (min_density < 0)
        nodes_kept = find(pi_orig > (1/n_nodes));
    else
        nodes_kept = find(pi_orig > min_density);
    end
    
else
    %TODO: add no sampling approach that auto determines epsilon
    disp('Sampling type not recognized.')
    return;
end
disp(['... done in ', num2str(toc), 's']);

n_nodes_kept = numel(nodes_kept);
data.NodeIndices = data.NonZeroIndices(nodes_kept);

W = W(nodes_kept, nodes_kept);

% normalize by power of density according to family parameter
d = full(sum(W,2));
inv_D_alpha = spdiags(d.^-p.alpha, 0, n_nodes_kept, n_nodes_kept);
W = inv_D_alpha*W*inv_D_alpha;



%% get eigenvectors of normalized laplacian
disp('Taking eigendecomposition of normalized laplacian ...'); tic;
[MsEigenvectors, MsEigenvalues] = normalized_laplacian_eigs(W, p.n_dims);
disp(['... done in ', num2str(toc), 's']);

%% writing values into map structure
map.Epsilon = p.epsilon;
map.CorrelationThreshold = p.correlation_thresh;
map.SignificanceThreshold = p.significance_thresh;
map.ThresholdType = p.thresh_type;
map.MinimumNodeDegree = p.min_corr_sum;
map.MinimumDensity = p.min_density;
map.SamplingType = p.sampling_type;
map.Alpha = p.alpha;

map.Degrees = d;
map.MsEigenvectors = MsEigenvectors;
map.MsEigenvalues = MsEigenvalues;