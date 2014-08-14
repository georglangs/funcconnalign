function [E V] = updatePhis_dd_danial(s,L,data,params,E,V,z,nos)

ordinds = randperm(nos.N(s));

degrees = data.degrees{s};
D = sparse((1:data.nos.N(s))',(1:data.nos.N(s))',degrees);

for k = 1:params.K
    muInvSig(:,k) = params.cluster.mu(:,k)'*params.cluster.invsig(:,:,k);
end
zMuInvSig = z * muInvSig';
for count = 1:nos.N(s)
    zInvSig(:,:,count) = sum( ...
        permute(z(count*ones(nos.D,1),:,ones(nos.D,1)),[1 3 2]) .* ...
        params.cluster.invsig ,3);
end

phiPhi = full(E'*D*E) + diag(V'*degrees);

for count = 1:nos.N(s)
    
    voxind = ordinds(count);
    
    phiPhi = phiPhi - degrees(voxind)*E(voxind,:)'*E(voxind,:);
    phiPhi = phiPhi - degrees(voxind)*diag(V(voxind,:));
    
    for d = 1:nos.D
        
        V(voxind,d) = params.rho(s)*degrees(voxind)*phiPhi(d,d) + ...
            z(voxind,:)*squeeze(params.cluster.invsig(d,d,:));
        V(voxind,d) = 1./V(voxind,d);
        
        %if V(voxind,d) < 0
        %    keyboard
        %end
        
        E(voxind,d) = degrees(voxind)*params.rho(s) ...
            * (    E(:,d)'*D*L.L(:,voxind)  ...
            -E(voxind,:)*phiPhi(:,d) + E(voxind,d)*phiPhi(d,d) ) ...
            +  zMuInvSig(voxind,d) -E(voxind,:)*zInvSig(:,d,voxind) ...
            + E(voxind,d)*zInvSig(d,d,voxind);
        E(voxind,d) = E(voxind,d)*V(voxind,d);
        
        
    end
    E = full(E);
    V = full(V);
    phiPhi = phiPhi + degrees(voxind)*E(voxind,:)'*E(voxind,:);
    phiPhi = phiPhi + degrees(voxind)*diag(V(voxind,:));
    
    1;
end
