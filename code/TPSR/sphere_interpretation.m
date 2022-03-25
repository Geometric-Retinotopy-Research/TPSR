%% 
% Description  --  function Vqdst = sphere_interpretation(Vsrc, Vdst, Vq)
%		interperating sphere.		A map between Vsrc to Vdst on sphere. 		
%       Then given query points on the sphere (Vsrc domain), we can
%       interpreate the Vq in the dst domain 
% Parameter(s): 
%     Vsrc[Nx3 double array]   --  Map source
%     Vdst[Nx3 double array]   --  Map dst
%     Vq[Mx3 double array]   --    Query point based on the Map
% Return:
% 	  Vqdst[Mx3 double array]  -- the result of interpretation
%%
function Vqdst = sphere_interpretation(Vsrc, Vdst, Vq)

 

% construct a map between vsrc to vdst, the compute the dest point for vq 
% assert(abs(norm(vsrc)-100))
idn = Vsrc(:,3)<1;   ids = Vsrc(:,3)> -1;

% north 
for i=1:3
    u = north_prj(Vsrc(idn,:));
    Fi = scatteredInterpolant(u, Vdst(idn,i));    
    uq = north_prj(Vq);    
    bd = boundary(u,0);
    in = inpolygon(uq(:,1), uq(:,2), u(bd,1), u(bd,2));    
    Vqdst(in,i) = Fi(uq(in,:));
end

% south
for i=1:3
    u = south_prj(Vsrc(ids,:));
    Fi = scatteredInterpolant(u, Vdst(ids,i));    
    uq = south_prj(Vq);
    bd = boundary(u,0);
    in = inpolygon(uq(:,1), uq(:,2),u(bd,1), u(bd,2));    
    Vqdst(in,i) = Fi(uq(in,:));
end



end


function u = north_prj(V)
u = [V(:,1)./(100-V(:,3)) V(:,2)./(100-V(:,3)) ];
end

function u = south_prj(V)
V(:,3) = -V(:,3);
u = [V(:,1)./(100-V(:,3)) V(:,2)./(100-V(:,3)) ];
end