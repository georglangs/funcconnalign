function [params, hvars, data, pointer] = perform_varEM_dd_fullcov(map, cluster, rho, file_name, params, hvars, data, pointer)
 
% In the new version, we subtract the projection of the L matrix into the subspace spanned by th
% constant eigenvector, so we deal with the first uninformative dime    nsion
% off-line.

if isfield(map,'diffusion_time')
diffusion_time = map.diffusion_time;
else % i.e. unknown
   diffusion_time = 0; 
end
rmpath mfiles
addpath mfilesDanial
addpath ../../matlab_packages/netlab.3.3
%addpath ../matlab_packages/VariationalBayesEM_GMM
addpath ../../matlab_packages/GeorgHelperFunctions/
%addpath /Volumes/BASIL/georg/CommonCode/matlab_packages/GeorgHelperFunctions

addpath /afs/csail.mit.edu/u/l/langs/georg/CommonCode/matlab_packages/GeorgHelperFunctions

% clear params hvars data pointer
%matlabpool

% Initialization:
%--------------------------------------------------------


mymap = map;

if ( nargin < 3)
    %mymap = jointmap;
    % we are using: map.Phi, map.L, map.clusterMu, map.clusterSigma
    data.nos.D = size(mymap.Phi,2)-1;   % dimension of embedding 
    data.nos.S = length(map.SubjectIndices); % number of subjects
    params.K = size(cluster.clusterPosterior,2);
    
     data.degrees_mean = zeros(data.nos.S,1);
    
    for counter = 1:data.nos.S
        data.L{counter}.L = mymap.L{counter};
        data.L{counter}.L = data.L{counter}.L - diag(diag(data.L{counter}.L));
        data.degrees{counter} = map.degrees{counter};
        data.degrees_mean(counter) = mean(data.degrees{counter}(:));                      % <---------------- d_i d_j
    end
    
    ss = 0;
    for s = 1:data.nos.S
        data.nos.N(s) = size(data.L{s}.L,1);
        %for count = 1:data.nos.N(s)
        %    data.L{s}.neighbs{count} = find(data.L{s}.L(count,:));
        %end
        data.L{s}.L = data.L{s}.L - mymap.Phi(mymap.SubjectIndices{s},1)*mymap.Phi(mymap.SubjectIndices{s},1)';
        data.L{s}.sum2 = sum(data.L{s}.L(:).^2);
        data.L{s}.Id = (data.L{s}.L > 0);
        pointer{s} = (ss+1:ss+data.nos.N(s)); % points to the points of one subject
        hvars.phi.E(pointer{s},:) = mymap.Phi(mymap.SubjectIndices{s},2:data.nos.D+1)  % Ntot X D
        ss = ss+ data.nos.N(s);
    end
    
    hvars.phi.V = 1e-4*ones(sum(data.nos.N),data.nos.D);   % Ntot X D
    %hvars.phi.V(:,1) = zeros(sum(data.nos.N),1);

    % Hidden Variables:

    %hvars.phi.E = mymap.Phi(:,1:data.nos.D)  % Ntot X D
    %hvars.phi.V = 1e-8*ones(sum(data.nos.N),data.nos.D);   % Ntot X D
    hvars.z = cluster.clusterPosterior;    % Ntot X K

    % Parameters

    params.cluster.mu = cluster.clusterMu(2:end,:);  % D X K
    params.pi = sum(cluster.clusterPosterior)/sum(data.nos.N);
    %{
    params.cluster.sig = cluster.clusterSigma(2:end,2:end,:); % D X D X K
    for count = 1:params.K
        params.cluster.invsig(:,:,count) = inv(cluster.clusterSigma(2:end,2:end,count));
    end
    params.pi = ones(params.K,1)';
    %}
    
    params = updateClusters(data,params,hvars);
    
%     params.rho = 10^7*ones(data.nos.S,1);
%     params.rho = 10^5*ones(data.nos.S,1);
%     params.rho = 10^5*ones(data.nos.S,1)/mean(data.degrees_mean);
    
    % params.rho = (rho ./ mean(cell2mat(map.degrees'))^2)*ones(data.nos.S, 1);
    params.rho = rho*ones(data.nos.S,1) ./ data.degrees_mean.^2;
    
    %for s=1:6,avdegrees=mean(atlas_nr.degrees{s});end,params.rho = 10%5./avdegrees;
    %params.rho = 10^8*ones(data.nos.S,1);
      
    
end

hvars0 = hvars;
params0 = params;
% Main Loop
%------------------------------------------------------------
hvars = hvars0;
params = params0;
%
condition = 1
counter = 1;
%
clear E
Gold = calculateGibbs_dd_danial(data,params,hvars,pointer);
%

ImprovementRatio = 10;

disp('here');
% save('init_temp','hvars','params','data')
%visuscript
%load('init_temp','hvars','params','data')

while (counter<31) %(abs(ImprovementRatio) >.000000001)|(counter<20)
% _dddd -> 0.0001
% _dddd_long -> 0.00001
% _dddd_long_ipmi -> 0.000001
% i have reduced it from .000001 -> .00001 for time reasons
% it was .0001 but georg set it to .00001 on miccai deadline day

    %clc
    disp(['Interation number = ' num2str(counter)]);
    %
    if counter > 2
        hvars = updateZs(data,params,hvars); % FIX
        disp('Z')
      %  G = calculateGibbs_dd_danial(data,params,hvars,pointer); sum(G-Gold)
      %  Gold = G;
    end
    
    %
    if counter > 10
        params = updateRho_dd_danial(data,params,hvars,pointer); % FIX            % <---------------- d_i d_j
        disp('rho')
      %  G = calculateGibbs_dd_danial(data,params,hvars,pointer); sum(G-Gold)
      %  Gold = G;
    end
    %}
    tic
    
    for s = 1:data.nos.S
     [hvars.phi.E(pointer{s},:) hvars.phi.V(pointer{s},:)] = ...  
            updatePhis_dd_danial(s,data.L{s},data,params, ...                  % <---------------- d_i d_j
            hvars.phi.E(pointer{s},:), hvars.phi.V(pointer{s},:),...
            hvars.z(pointer{s},:),data.nos);
    end
    toc
    disp('phi')
    %G = calculateGibbs_dd_danial(data,params,hvars,pointer); sum(G-Gold)
    %Gold = G;
    
    tic;
    params = updateClusters(data,params,hvars);
    toc
    disp('clusters')
    G = calculateGibbs_dd_danial(data,params,hvars,pointer);
    % sum(G-Gold)
    %Gold = G;
    
    hvars.GibbsFreeEnergy(counter,:) = G;
    
    if counter>3
    ImprovementRatio = sum(hvars.GibbsFreeEnergy(counter-1,:)-hvars.GibbsFreeEnergy(counter,:),2) / sum(hvars.GibbsFreeEnergy(counter-1,:));
    else
        ImprovementRatio = 1;
    end
   if mod(counter,5)==0
%        save(['IntermedResult_VarEM_debugdanial_t' num2str(diffusion_time) '_K' num2str(params.K) '_D' num2str(data.nos.D) 'It' num2str(counter)], ...
%            'hvars','params','pointer');
        save([file_name, '_iter' num2str(counter), '.mat'], '-mat', 'hvars','params','pointer');
   end
    %params.pi = updatePi(data,params,hvars)
    %visuscript
    
    
    mean(hvars.phi.V)
    
    %  visuscript, title('pi') ; drawnow; %pause
    counter = counter +1
end
%  save(['IntermedResult_VarEM_t' num2str(diffusion_time) '_K' num2str(params.K) '_D' num2str(data.nos.D) 'It' num2str(counter) '_Final'], ...
%            'hvars','params','pointer','hvars0','params0'); 