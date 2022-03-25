 %%
% Description  -- function data =read_nii_info(fn)
%       read .nii file
% 
% Parameter(s):
% 		fn[string] -- file name.
% return:
%       data[double array]  -- content of the file
%
%%
% Smooth on the boundary  by line fit
function bd_pos = smooth_bd_vis_byfit(uv_pol, vis,  bd, R2)

rho = uv_pol(:,1);
theta = uv_pol(:,2);
bd_pos = zeros(length(bd),2);
R2(isnan(R2))=0;
for i=1:length(bd)     % for each bd
    
    
    bdi  = bd(i);
    featurei = uv_pol(bdi,1); % i-th boundary rho
    
    % find points in close to line roh
    th_fd = 0.1;
    while 1
        id = find(abs(rho - featurei)<th_fd);
        if length(id) <10
            th_fd = th_fd*1.2;
        else
            break;
        end
    end
    
    [thetas, sid ]= sort(theta(id));
    sorted_id = id(sid);
    
    
    value_i = vis(sorted_id,2);
    
    % estimate Weigth R2line
    R2line = 20*weight_from_R2(R2(sorted_id));
    R2line(abs(diff(value_i))>0.5)=0;
    mm = movmean(R2line,8);
    R2line(mm<0.6)=0;
    
    
    % smooth ang
    value_smoothed = smoothn(value_i, R2line+0.01);
    s = spline(thetas,value_smoothed,theta(bdi));
    
    bd_pos(i,2) =s;
end



for i=1:length(bd)
    
    bdi  = bd(i);
    featurei = uv_pol(bdi,2); % i-th boundary rho
    
    % find points in close to line roh
    
    th_fd = 0.1;
    while 1
        id = find(abs(theta - featurei)<th_fd);
        if length(id) <10
            th_fd = th_fd*1.2;
        else
            break;
        end
    end
    
    [rhos, sid ]= sort(rho(id));
    sorted_id = id(sid);
    
    
    value_i = vis(sorted_id,1);
    
    % estimate Weigth R2line
    R2line = 20*weight_from_R2(R2(sorted_id));
    R2line(abs(diff(value_i))>0.5)=0;
    mm = movmean(R2line,8);
    R2line(mm<0.6)=0;
    
    
    % smooth ecc
    value_smoothed = smoothn(value_i, R2line+0.01);
    
    tx = rhos;
    ty = value_smoothed;
    
    [tx, ia] = unique(tx);
    ty = ty(ia);
    s = spline(tx, ty, rho(bdi));    
    bd_pos(i,1) =s;
    
end


x = bd_pos(:,1);
xe = smoothn([x;x;x],[R2(bd);R2(bd);R2(bd)]);
bd_pos(:,1) = xe(length(x)+1:2*length(x));


x = bd_pos(:,2);
xe = smoothn([x;x;x],[R2(bd);R2(bd);R2(bd)]);
bd_pos(:,2) = xe(length(x)+1:2*length(x));
% bd_pos= prune_bd_vis(bd_pos);

end