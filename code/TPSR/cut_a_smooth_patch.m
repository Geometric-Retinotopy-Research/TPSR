%%
% Description  -- function smoothed = curve_smooth(anchorpos, R2)
%		cut a smooth patch start from a boundary described by polygon
% Parameter(s):
%       Froi[double array]	        --  Face
%       Vroi[String]                --  verticles
%       polygon                     --  the amount of sides
% return:
%		smoothed[double array]    --  results of smoothed
%
%%
% 
function [Fcut, Vcut, vfathercut]= cut_a_smooth_patch(Froi, Vroi, uvroi, polygon)
[in, on] = inpolygon(uvroi(:,1),uvroi(:,2),polygon(:,1),polygon(:,2));
Nv = size(uvroi,1);
id2delete = setdiff(1:Nv, find(in|on));
[Fcut, Vcut, vfathercut] = gf_remove_mesh_vertices(Froi, Vroi, id2delete);

 
% 
% bd = compute_bd(Fcut);
% bdpos = uvroi(vfathercut(bd),:);
% 
% bdpos_smooth = [ext_smooth(bdpos(:,1)) ext_smooth(bdpos(:,2))];
% dis = sqrt(dot(bdpos_smooth - bdpos_smooth,bdpos_smooth - bdpos_smooth, 2));
% fixids = find( dis< 0.3*mean(dis));
% 
% bdpos_smooth(fixids,:) = bdpos(fixids,:);
% 
% [in, on] = inpolygon(uvroi(:,1),uvroi(:,2),bdpos_smooth(:,1),bdpos_smooth(:,2));
% Nv = size(uvroi,1);
% id2delete = setdiff(1:Nv, find(in|on));
% [Fcut, Vcut, vfathercut] = gf_remove_mesh_vertices(Froi, Vroi, id2delete);
%  figure;plot_mesh(Fcut,uvroi(vfathercut,:))
%  
 
end

function sx = ext_smooth(x)
x =x(:);
n =length(x);
sxe = smooth([x;x;x],3);
sx = sxe(n+1:2*n);
end