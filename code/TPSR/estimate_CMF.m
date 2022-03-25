%%
% Description  -- function CM = estimate_CMF(F,V,visxy, k)
%					return a estimate CMF value;
% Parameter(s):
%       F[double array]	  --  face 
%       V[douhle array]   --  vertices
%		vis[double array] --  visual area
% 
% return:
%		CM[double array]  --  return a CMF value
%
%%
function CM = estimate_CMF(F,V,visxy, k)
if nargin<4
    k=2;
end

vrs =compute_vertex_kring(F,V,1:size(V,1), k);

va = vertex_area(F,V);

% cortical area
n = size(visxy,1);
CM = zeros( n,1);
for i=1:n
    cortical_area = sum(va(vrs{i}));
    
    rth =visxy(vrs{i},:);  % be aware it is polar system
    
    [~,visual_area] = boundary([rth(:,1).*cos(rth(:,2)) rth(:,1).*sin(rth(:,2))]);     
    
    CM(i) = cortical_area/visual_area;
end

end
