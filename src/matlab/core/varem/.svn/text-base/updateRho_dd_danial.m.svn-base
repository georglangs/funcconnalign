function params = updateRho_dd_danial(data,params,hvars,pointer)


	for s = 1:data.nos.S

	    phi = hvars.phi.E(pointer{s},:);
	    V = hvars.phi.V(pointer{s},:);
	    phi2 = phi.^2;
		D = sparse((1:data.nos.N(s))',(1:data.nos.N(s))',data.degrees{s});
		L = data.L{s}.L;

	    params.rho(s) = ...
	    data.degrees{s}'*(L.^2)*data.degrees{s} - 2*trace(phi'*D*L*D*phi)+trace(phi'*D*phi*phi'*D*phi) ...
	               - sum((sum(phi2').*data.degrees{s}').^2) ...
	               +2*( data.degrees{s}'*phi2*V'*data.degrees{s} -trace(D*phi2*V'*D) ) ...
	               + data.degrees{s}'*V*V'*data.degrees{s} - sum(sum(V.^2,2).*data.degrees{s});

	               %{
	               	rhoval(s) = ...
    trace(D*L*D*L) - 2*trace(phi'*D*L*D*phi)+trace(phi'*D*phi*phi'*D*phi) ...
               - sum((sum(phi2').*data.degrees{s}').^2) ...
               +2*( data.degrees{s}'*phi2*V'*data.degrees{s} -trace(D*phi2*V'*D) ) ...
               + data.degrees{s}'*V*V'*data.degrees{s} - sum(sum(V.^2,2).*data.degrees{s});
	               	%}
	               
	               
	    params.rho(s) = data.nos.N(s)*(data.nos.N(s) -1)/full(params.rho(s));

	end


