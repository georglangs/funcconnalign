function hvars = updateZs(data,params,hvars)

for k = 1:params.K
    
    logz(:,k) = -.5*sum((hvars.phi.E*params.cluster.invsig(:,:,k)).*hvars.phi.E,2) ...
                -.5*params.cluster.mu(:,k*ones(sum(data.nos.N),1))'*params.cluster.invsig(:,:,k)*params.cluster.mu(:,k) ...
                -.5*hvars.phi.V*diag(params.cluster.invsig(:,:,k)) ...
                + hvars.phi.E*params.cluster.invsig(:,:,k)*params.cluster.mu(:,k) ...
                -.5*ones(sum(data.nos.N),1)*(log(det(params.cluster.sig(:,:,k))) + logbz(params.pi(k)));
    
end
logz = logz - max(logz,[],2)*ones(1,params.K);
z = exp(logz);
sumz = sum(z,2);
hvars.z = z ./ (sumz*ones(1,params.K));
