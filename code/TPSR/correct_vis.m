%%
% Description  -- function vis = correct_vis(Em, lr)
%		return the result if correct with vis EM and lr
% Parameter(s):
%       Em[double array]	  --  
%       lr[String]            --  type of the hemisphere.
% return:
%		vis[double array]     --  return visual area.
%
%%
function vis = correct_vis(Em, lr)
prf_copy = Em.Vertex_prf;
if nargin < 2
%     disp('warning: no lr We assume left hemisphere');
     lr='lh';
end
if (strcmp(lr, 'lh'))  
%   Columns 1 through 11
% 
%     {'Unknown'}    {'V1v'}    {'V1d'}    {'V2v'}    {'V2d'}    {'V3v'}    {'V3d'}    {'hV4'}    {'VO1'}    {'VO2'}    {'PHC1'}
% 
%   Columns 12 through 22
% 
%     {'PHC2'}    {'TO2'}    {'TO1'}    {'LO2'}    {'LO1'}    {'V3B'}    {'V3A'}    {'IPS0'}    {'IPS1'}    {'IPS2'}    {'IPS3'}
% 
%   Columns 23 through 26
% 
%     {'IPS4'}    {'IPS5'}    {'SPL1'}    {'FEF'}
    prf_copy(:,1) = prf_copy(:,1)/180*pi; 
    prf_copy(: ,1) = mod(prf_copy(: ,1) + pi, 2*pi );
    
    prf_copy(Em.Vertex_atlaswang==3 ,1) = -prf_copy(Em.Vertex_atlaswang==3 ,1) + 3*pi;
    prf_copy(Em.Vertex_atlaswang==4 ,1) = -prf_copy(Em.Vertex_atlaswang==4 ,1) + pi;
%     
    prf_copy(Em.Vertex_atlaswang==5 ,1) =  prf_copy(Em.Vertex_atlaswang==5 ,1) + pi;
    prf_copy(Em.Vertex_atlaswang==6 ,1) =  prf_copy(Em.Vertex_atlaswang==6 ,1) - pi;    
    
    prf_copy(Em.Vertex_atlaswang==16 ,1) =  prf_copy(Em.Vertex_atlaswang==16 ,1) +  2*pi;  
else
    
    error('To Be Done');
    
end
 
vis = prf_copy(:,[2,1]); % in the order of r, theta

end