

function [hat_vis, flip] = topological_smoothing_test3(face, uv_pol, vis, R2, bd_id, bd_vis0, bd_tol, smooth_lambda0, smooth_avg_k, meanddth)

if nargin < 10
    meanddth = 1.0;
end
if nargin <9        
    smooth_avg_k = 2;
end

if nargin <8    
    smooth_lambda0 = 0.001;
end

if nargin <7    
    bd_tol = 0;
end
 

%% Average Smoother
  
vr = compute_vertex_ring(face);

fixed  =0 ;
for threshold = 1:-0.02:0.05
    hat_vis = vis;
    W = R2;
    W(bd_id) = W(bd_id)+10;

    for si=1:20        
        [flip, fmu, vmu] = get_flips_nobd(face, uv_pol, hat_vis);
        nonsmoothid = find(abs(vmu)>threshold)';
        hat_vis = avg_smooth(hat_vis, vr,nonsmoothid , W);
        W(nonsmoothid) = W(nonsmoothid)*0.9;
        W(bd_id) = W(bd_id)*1.111;
%         fprintf('Threshold = %f Iteration %d : Flip %d\n', threshold, si, flip)
        if flip ==0
            fixed = 1;
            break
        end
    end
    meanse = nanmean(vecnorm( (hat_vis-vis)'));
    if meanse>0.5      
        continue
    end
    if fixed
        break
    end
end


end

function vis = avg_smooth(vis, vr, roi, W)

    vis_avg = vis;
    for i = roi
        vri = vr{i};
        vis_avg(i,:) =  sum([vis(vri,1).*W(vri)  vis(vri,2).*W(vri)],1)/sum(W(vri));
    end
    vis = vis_avg;
end

 
% compute the number of flips without bd
function [flip, fmu, vmu] = get_flips_nobd(face, uv, fn)
% bd = compute_bd(face);
% fv = face_area(face,fn); 
% fu = face_area(face,uv); 
% face_sim = face;
% face_sim(ismember(face(:,1),bd)| ismember(face(:,2),bd) |ismember(face(:,3),bd) | fv<1e-5| fu<1e-5,: )=[]; 
flip = length(find(abs(compute_bc(face,uv,fn))>1));
fmu = compute_bc(face,uv,fn);
fmu(abs(fmu)>0.99) = fmu(abs(fmu)>0.99) * 20;
Operator = createOperator(face,uv);
vmu = Operator.f2v*fmu; 
end




%% rect bd_vis
function bd_vis_r = rect_vis(bd_vis, bd_vis0, bd_tol)
bd_vis_r = bd_vis;
for i=1:size(bd_vis,1)
    ni = norm(bd_vis(i,:) - bd_vis0(i,:));
    if  ni> bd_tol
        bd_vis_r(i,:) = bd_vis0(i,:) + (bd_vis0(i,:) - bd_vis(i,:))/ni * bd_tol;
    end
end

end