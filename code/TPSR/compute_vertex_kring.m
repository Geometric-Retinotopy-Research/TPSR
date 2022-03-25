%%
% Description  -- function vkr = compute_vertex_kring(face,vertex,vc,K, rotate)
%
% Parameter(s):
%       face[double array]    --  connectivity of mesh.
%       vertex[double array]  --  vertex of mesh.
%		vc[double array]      --  
%		K[double array]       -- 
%       rotate                -- rotate information.
% return:
%		vkr[double array]     -- vertex kring.
%
%%
function vkr = compute_vertex_kring(face,vertex,vc,K, rotate)


if ~exist('rotate','var')  
    rotate = 0;
end

vr = compute_vertex_ring(face,vertex,[],1);

for i = 1: length(vc)
     
    neigh1=vr{vc(i)};
    
    rotate = mod(rotate,length(neigh1));
    neigh1 = neigh1([rotate+1:end 1:rotate]);
    
    PatchPts=[vc(i)  neigh1];
     
    
    for k = 2: K
        
        Pts2Add = [];
        for l=1:length(PatchPts)
               
           Pts2Add =[Pts2Add  vr{PatchPts(l)}];
        
        end
        
        PatchPts = [PatchPts  unique(Pts2Add,'stable')];   
        
        PatchPts = unique(PatchPts,'stable');
    
    end
   
    vkr{i} =  unique(PatchPts,'stable');
    
end


end