function params = updateClusters(data,params,hvars)

params.pi = sum(hvars.z);
params.cluster.mu = hvars.phi.E' * hvars.z ./ (ones(data.nos.D,1)*params.pi);

for k = 1:params.K
    params.cluster.sig(:,:,k) = ( hvars.z(:,k*ones(data.nos.D,1)).*( hvars.phi.E-params.cluster.mu(:,k*ones(sum(data.nos.N),1))' ) )' ...
                       * ( hvars.phi.E-params.cluster.mu(:,k*ones(sum(data.nos.N),1))' ) + diag(hvars.z(:,k)'*hvars.phi.V);
    params.cluster.sig(:,:,k) = params.cluster.sig(:,:,k) / params.pi(k);
    params.cluster.invsig(:,:,k) = inv(params.cluster.sig(:,:,k) );
end

params.pi =  params.pi / sum(params.pi); 


