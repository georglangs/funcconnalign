function Fval = calculateGibbs_dd_danial(data,params,hvars,pointer)
%% % Likelihood of affinities


for s = 1:data.nos.S
    s
    phi = hvars.phi.E(pointer{s},:);
    V = hvars.phi.V(pointer{s},:);
	D = sparse((1:data.nos.N(s))',(1:data.nos.N(s))',data.degrees{s});
	L = data.L{s}.L;
    phi2 = phi.^2;
    
    rhoval(s) = ...
    data.degrees{s}'*(L.^2)*data.degrees{s} - 2*trace(phi'*D*L*D*phi)+trace(phi'*D*phi*phi'*D*phi) ...
               - sum((sum(phi2').*data.degrees{s}').^2) ...
               +2*( data.degrees{s}'*phi2*V'*data.degrees{s} -trace(D*phi2*V'*D) ) ...
               + data.degrees{s}'*V*V'*data.degrees{s} - sum(sum(V.^2,2).*data.degrees{s});
				%before:
				%rhoval(s) = ...
				%trace(D*L*D*L) - 2*trace(phi'*D*L*D*phi)+trace(phi'*D*phi*phi'*D*phi) ...
				%- sum((sum(phi2').*data.degrees{s}).^2) ...
				%+2*( data.degrees{s}'*phi2*V'*data.degrees{s} -trace(D*phi2'*V*D ) ...
				%+ data.degrees{s}'*V*V'*data.degrees{s} - sum(sum(V.^2,2).*data.degrees{s}); ...
																	
	rhoval(s) = full(rhoval(s));           
    rhoval(s) = 0.25*params.rho(s)*rhoval(s)-0.25*data.nos.N(s)*(data.nos.N(s) -1)*log(params.rho(s));
    
end

% Likelihood of affinities
Fval(1) = sum(rhoval);

%% Lieklihood of embdedded vectors

Fval(2) = 0;
for k = 1:params.K
    sigval = ( hvars.z(:,k*ones(data.nos.D,1)).*( hvars.phi.E-params.cluster.mu(:,k*ones(sum(data.nos.N),1))' ) )' ...
                       * ( hvars.phi.E-params.cluster.mu(:,k*ones(sum(data.nos.N),1))' ) + diag(hvars.z(:,k)'*hvars.phi.V);
    Fval(2) = Fval(2) + 0.5*sum(sum(sigval.*params.cluster.invsig(:,:,k))) + 0.5*sum(hvars.z(:,k))*logbz(det(params.cluster.sig(:,:,k)));
end	

%% prior on zs

Fval(3) = -sum(hvars.z)*logbz(params.pi)';

%% Entropy of z:

Fval(4) = hvars.z(:)'*logbz(hvars.z(:));

%% Entropy of phi

Fval(5) = -0.5* sum(logbz(hvars.phi.V(:)));

