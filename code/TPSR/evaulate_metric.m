%%
% Description  -- function [meanse, stdse, mean_ang, std_ang, flip]=evaulate_metric(fn,f,face,uv)
%		 evaulate results of smoothing.
%
% Parameter(s):
%		subjectfn[String]     --  file name, read data from this file.
%		lr[String]            --  type of hemisphere.
%		hat_vis[double array] --  x,y in pixel.
%		uvw[double array]     --  uv coordinary
% return:
%		NRMSE_raw[double]   --  Normalized RMSE of the pRF solution in fMRI	fitting.
%		NRMSE_new[double]   --  Normalized RMSE of the pRF0 solution in fMRI fitting
%		R2_raw[double] --  Average Valiance explained for pRF solution
%		R2_new[double] --  Average Valiance explained for pRF0 solution
%		AIC_raw[double]     --  AIC metric for pRF solution
%		AIC_new[double]     --  AIC metric for pRF0 solution
%		extra[double]       --  extra information.
%
%%
function [meanse, stdse, mean_ang, std_ang, flip]=evaulate_metric(fn,f,face,uv)
global ang
meanse = nanmean(vecnorm( (fn-f)'));
stdse = nanstd(vecnorm( (fn-f)'));

ang = estimate_angle(face,uv,fn,linspace(min(fn(:,1)), max(fn(:,1)), 10) , linspace(min(fn(:,2)), max(fn(:,2)), 10) );

ang(isnan(ang))=[];

mean_ang =nanmean(abs(ang-pi/2));
std_ang = nanstd(abs(ang-pi/2));



% % remove face lay on boundary
bd = compute_bd(face);
fv = face_area(face,fn); 
fu = face_area(face,uv); 
face(ismember(face(:,1),bd)| ismember(face(:,2),bd) |ismember(face(:,3),bd) | fv<1e-4| fu<1e-4,: )=[]; 
flip = length(find(abs(compute_bc(face,uv,fn))>1)); 

end