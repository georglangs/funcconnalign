function [map,W] = atlas_build_functional_map(data, p, verbose)
% Function to build a functional geometry map from fMRI data
% *** This is a preliminary release, please don't distribute ***
%
% INPUT:
% data.BOLD ... BOLD signal data, each row corresponds to one voxel, each
% column to one time point
% optional alternative to BOLD:
% data.affinitymatrix ... if precomputed, then you can provide the affinity
% matrix directly. no need for BOLD in this case.
%
% Georg Langs
% 23.9.2011
%

if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end

%%
if isfield(data,'affinity_matrix')
  p.do_only_embedding = 1;
else
  p.do_only_embedding = 0;
end

%%
if p.do_only_embedding == 0
data.NonZeroIndices = max(abs(data.BOLD),[],2)>0;
data.NonZeroIndicesNumbers = find(data.NonZeroIndices);
else
data.NonZeroIndices = ones(size(data.affinity_matrix,1),1);
data.NonZeroIndicesNumbers = find(data.NonZeroIndices);
end
%% Set parameters:


dp.do_save_memory = 1;
% 0 ... create corr coef matrix densely then create sparse
% 1 ... create corr coef matrix sparsely (saves memory at expense of proc time)
do_define_gamma = 1;
% 0 ... nothing
% 1 ... define diffusion map coordinates

dp.map_dimension = 20;
dp.delta = [];
dp.diffusion_time = 2;

dp.n_eigs = 50;
dp.epsilon = -1;

dp.min_corr_sum = 100;
dp.min_density = -1; % if negative then automatically choose better than chance

dp.keep_negative_correlations = 0;
dp.alpha = 0;


dp.thresh_type = 'significance';
dp.significance_thresh = 0.01;

% dp.thresh_type = 'correlation';
% dp.correlation_thresh = 0.6;

% dp.sampling_type = 'corr_sum'; 
dp.sampling_type = 'all';

% parameters during session with satra 13.12.2013:
dp.do_exp_rescaling_of_correlations = 1;
dp.do_alpha_normalization = 1;
dp.do_diffusion_operator_eigs = 0;


p = defaultParams(p,dp);




%% Create correlation coefficient matrix
if p.do_only_embedding == 0
    if verbose
        disp('Generating correlation matrix ...');
    end
    tic;
    if p.do_save_memory
        
        if (strcmp(p.thresh_type, 'significance'))
            GCT = atlas_sparse_corrcoef(data.BOLD(data.NonZeroIndices,:)', p.significance_thresh, p.thresh_type, p.keep_negative_correlations, verbose);
            
        elseif (strcmp(p.thresh_type, 'correlation'))
            GCT = atlas_sparse_corrcoef(data.BOLD(data.NonZeroIndices,:)', p.correlation_thresh, p.thresh_type, p.keep_negative_correlations, verbose);
        else
            error('No other threshold types available yet. ')
        end
        
    else
        GC = corrcoef(data.BOLD(data.NonZeroIndices,:)');
        
        if (strcmp(p.thresh_type, 'correlation'))
            GCT = sparse(GC.*(GC > p.correlation_thresh));
        else
            error('No other threshold types available yet. ')
        end
        
        GCT(isnan(GCT)) = 0;
        clear GC
    end
    if verbose
        disp(['... done in ', num2str(toc), 's']);
    end
else
    tic;
    if verbose
        disp('Using pre computed correltaion matrix, performing sparcification ...');
    end
    GC = data.affinity_matrix;
    if (strcmp(p.thresh_type, 'correlation'))
        GCT = sparse(GC.*(GC > p.correlation_thresh));
    else
        error('No other threshold types available yet. ')
    end
    
    GCT(isnan(GCT)) = 0;
    clear GC
    if verbose
        disp(['... done in ', num2str(toc), 's']);
    end
end

%% Calculate the degree of the nodes, and discard nodes with degree below threshold
if verbose
    disp('Generating and sampling affinity matrix ...'); tic;
end
n_nodes = numel(data.NonZeroIndices);

% automatically determine epsilon if required
if ( p.epsilon < 0 )
    epsilon_used = median(nonzeros(GCT));
    if verbose
        disp(['Determined epsilon = ', num2str(epsilon_used), ' , from median.']);
    end
else
    epsilon_used = p.epsilon;
end

if verbose                           
    disp('Creating W ...'); tic;
end
%W = sparse(n_nodes, n_nodes);
%if verbose
%    disp('... initialized, next exp() ...'); 
%end
if p.do_exp_rescaling_of_correlations == 1
 W = spfun(@(x) x/epsilon_used,GCT);
 W = spfun(@exp,W);
else
    disp('No rescaling performed i.e, W = GCT, epsilon is obsolete (atlas_build_functional_map line 152)')
   W = GCT; 
end 
%W(GCT~=0) = exp(GCT(GCT~=0) / epsilon_used);
if verbose
    disp('... done, now calculatin degrees	 ...'); 
end
graph_d = full(sum(GCT~=0,1)/2)';
if verbose
    disp(['... done in ', num2str(toc), 's']);
end
if verbose
    disp('Sub-sampling W ...'); tic;
end
if (strcmp(p.sampling_type, 'corr_sum'))
    
    corr_sums = full(sum(GCT));
    nodes_kept = find(corr_sums > min_corr_sum);
    
elseif (strcmp(p.sampling_type, 'density'))
    
    d_orig = full(sum(W,2));
    pi_orig = d_orig / sum(d_orig);
    
    if (p.min_density < 0)
        nodes_kept = find(pi_orig > (1/n_nodes));
    else
        nodes_kept = find(pi_orig > min_density);
    end
elseif  (strcmp(p.sampling_type, 'all'))
    nodes_kept = 1:size(GCT,1);
    disp('... keeping all nodes')
else
    %TODO: add no sampling approach that auto determines epsilon
    nodes_kept = 1:size(GCT,1);
    error('Sampling type not recognized.')
end
if verbose
    disp(['... done in ', num2str(toc), 's']);
end

n_nodes_kept = numel(nodes_kept);
map.NodeIndices = data.NonZeroIndices(nodes_kept);
map.indices_of_input_node_encoded_in_the_map = data.NonZeroIndicesNumbers(nodes_kept);
W = W(nodes_kept, nodes_kept);

if p.do_alpha_normalization == 1
% normalize by power of density according to family parameter
d = full(sum(W,2));
inv_D_alpha = spdiags(d.^-p.alpha, 0, n_nodes_kept, n_nodes_kept);
W = inv_D_alpha*W*inv_D_alpha;
end
%% Get eigenvectors of normalized laplacian
if verbose
    disp('Taking eigendecomposition of normalized laplacian ...'); tic;
end

if p.do_diffusion_operator_eigs == 1
[MsEigenvectors, MsEigenvalues] = atlas_diffusion_operator_eigs(W, p.n_eigs);
else
[MsEigenvectors, MsEigenvalues] = atlas_normalized_laplacian_eigs(W, p.n_eigs);
end

if verbose
    disp(['... done in ', num2str(toc), 's']);
end

%% Write values into map structure
map.Epsilon = epsilon_used;

if (strcmp(p.thresh_type, 'significance'))
    map.SignificanceThreshold = p.significance_thresh;
elseif (strcmp(p.thresh_type, 'correlation'))
    map.CorrelationThreshold = p.correlation_thresh;
end

map.ThresholdType = p.thresh_type;
map.MinimumNodeDegree = p.min_corr_sum;
map.MinimumDensity = p.min_density;
map.SamplingType = p.sampling_type;
map.Alpha = p.alpha;
map.GraphDegrees = graph_d;
map.Degrees = d;
map.MsEigenvectors = MsEigenvectors;
map.MsEigenvalues = MsEigenvalues;


%% 4.) Create Functional Maps

[map.Gamma, determined_param] = atlas_generate_diffusion_map(map.MsEigenvectors, map.MsEigenvalues, map.Degrees, p.map_dimension, p.diffusion_time, p.delta);

if (isempty(p.diffusion_time))
    if verbose
        disp(['Determined diffusion time  = ', num2str(determined_param)]);
    end
    map.DeterminedDiffusionTime = determined_param;
elseif (isempty(p.map_dimension))
    if verbose
        disp(['Determined map dimension  = ', num2str(determined_param)]);
    end
    map.DeterminedMapDimension = determined_param;
elseif (isempty(p.delta))
    if verbose
        disp(['Determined delta  = ', num2str(determined_param)]);
    end
    map.DeterminedDelta = determined_param;
end

map.Gamma(1,1)

map.Delta = p.delta;
map.DiffusionTime = p.diffusion_time;
map.MapDimension = p.map_dimension;
