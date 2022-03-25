
function [hat_vis, flip] = topological_smoothing_test2(face, uv_pol, vis, R2, bd_id, bd_vis0, changetol, smooth_lambda0, smooth_avg_k, meanddth)
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
    changetol = 0;
end
%% initilize hat_v 
hat_vis = vis;    

%%  
W = weight_from_R2(R2);
hat_vis = laplacian_smooth(uv_pol, hat_vis, W);

topo_flag = check_toplogical(face, uv_pol, hat_vis);

figure;
plot_mesh(face, hat_vis, double(topo_flag.* (R2>40)) );


end


%%
function topo_flag = check_toplogical(face, uv_pol, hat_vis)
vr = compute_vertex_face_ring(face, uv_pol);
muf =  compute_bc(face, uv_pol, hat_vis);
nv = size(uv_pol,1);
topo_flag =  zeros(nv,1);
for i=1:length(vr)
    vri = vr{i};
    if all(abs(muf(vri))<0.95)
       topo_flag(i)  = 1;
    end    
end
end

%% make_diffeomorphic:  try to make diffeomorphic by smooth and chop mu
function [hat_vis, flip] = make_diffeomorphic(face, uv_pol, hat_vis, bd_id, bd_vis)
oldhat = hat_vis;
% set parameters and operator
P.angle = 0.1;
P.upperBound = 0.999;
P.chopVal = 0.999;
P.smooth = 0.1;
Operator = createOperator(face,uv_pol);

% Smoothing and Chop on mu
for i = 1:5
    % compute mu
    mu = compute_bc(face,uv_pol, hat_vis);
    mu(isnan(mu))=0;
    
    % Smooth and chop on mu
    vmu = Operator.f2v*mu;
    Smooth_Operator = P.angle*speye(length(uv_pol)) + ...
        speye(length(uv_pol)) - P.smooth* Operator.laplacian/2;
    
    nvmu = Smooth_Operator\abs(vmu);
    vmu = nvmu.*(cos(angle(vmu))+sqrt(-1)*sin(angle(vmu)));
    smooth_mu = Operator.v2f*vmu;
    choped_mu = mu_chop(smooth_mu,P.upperBound,P.chopVal);
    
    % recovery function according to the smoothed and chopped mu
    [lbs_vis,~] = linear_beltrami_solver(face,uv_pol,choped_mu, bd_id, bd_vis);
    
    hat_vis = lbs_vis;
    
    % reverse back if nan occurs
    id = [find(isnan(hat_vis(:,1))); find(isnan(hat_vis(:,2)))];
    hat_vis(id,:) = oldhat(id,:) ;
    
    flip = get_flips(face, uv_pol, hat_vis);
    if( flip < 1 )
        break  % break when diffeomorphic
    end
end

end


function [meanvd, mean_dd] = visual_deviation(hat_vis, vis)

vis_distance =sqrt( (hat_vis(:,1).* (hat_vis(:,2) - vis(:,2) )).^2 + (hat_vis(:,1) - vis(:,1) ).^2);
meanvd = mean(vis_distance);
mean_dd = mean(vecnorm( (hat_vis-vis)'));
end

% compute the number of flips without bd
function flip = get_flips_nobd(face, uv, fn)
bd = compute_bd(face);
face(ismember(face(:,1),bd)| ismember(face(:,2),bd) |ismember(face(:,3),bd),:)=[]; 
flip = length(find(abs(compute_bc(face,uv,fn))>1));
end


% compute the number of flips
function [flip, flipid] = get_flips(face, uv, fn)
flipid = find(abs(compute_bc(face,uv,fn))>1);
flip = length(flipid);
end

% Smooth on the boundary locally
% kspan is the average over the span, take bigger to get more smoothed 3,4,5,...
function bd_pos = smooth_bd_vis_bybd(bd, bd_vis, kspan, flipidv)

n =size(bd_vis,1);
bd_visext = [bd_vis;bd_vis;bd_vis];
for i = 1:n
    bdi = bd(i);
    if ismember(bdi, flipidv)
        % if flip related boundary smooth locally
        bd_visext(n + i-5:i+ 5 + n,1) = smooth(bd_visext(n + i-5:i+ 5 + n,1), kspan);
        bd_visext(n + i-5:i+ 5 + n,2) = smooth(bd_visext(n + i-5:i+ 5 + n,2), kspan);
        
    else
        continue;
    end
end

bd_pos = bd_visext(n + 1:2*n,:);

end




% % Smooth on the boundary  by line fit
% function bd_pos = smooth_bd_vis_byfit( uv_pol, vis,  bd, R2, bd_vis0, changetol)
% 
% rho = uv_pol(:,1);
% theta = uv_pol(:,2);
% bd_pos = zeros(length(bd),2);
% 
% for i=1:length(bd)     % for each bd
%     
%     
%     bdi  = bd(i);
%     featurei = uv_pol(bdi,1); % i-th boundary rho
%     
%     % find points in close to line roh
%     th_fd = 0.1;
%     while 1
%         id = find(abs(rho - featurei)<th_fd);
%         if length(id) <10
%             th_fd = th_fd*1.2;
%         else
%             break;
%         end
%     end
%     
%     [thetas, sid ]= sort(theta(id));
%     sorted_id = id(sid);
%     
%     
%     value_i = vis(sorted_id,2);
%     
%     % estimate Weigth R2line
%     R2line = 20*weight_from_R2(R2(sorted_id));
%     R2line(abs(diff(value_i))>0.5)=0;
%     mm = movmean(R2line,8);
%     R2line(mm<0.6)=0;
%     
%     
%     % smooth ang
%     value_smoothed = smoothn(value_i, R2line+0.01);
%     s = spline(thetas,value_smoothed,theta(bdi));
%     
%     bd_pos(i,2) =s;
% end
% 
% 
% 
% for i=1:length(bd)
%     
%     bdi  = bd(i);
%     featurei = uv_pol(bdi,2); % i-th boundary rho
%     
%     % find points in close to line roh
%     
%     th_fd = 0.1;
%     while 1
%         id = find(abs(theta - featurei)<th_fd);
%         if length(id) <10
%             th_fd = th_fd*1.2;
%         else
%             break;
%         end
%     end
%     
%     [rhos, sid ]= sort(rho(id));
%     sorted_id = id(sid);
%     
%     
%     value_i = vis(sorted_id,1);
%     
%     % estimate Weigth R2line
%     R2line = 20*weight_from_R2(R2(sorted_id));
%     R2line(abs(diff(value_i))>0.5)=0;
%     mm = movmean(R2line,8);
%     R2line(mm<0.6)=0;
%     
%     
%     % smooth ecc
%     value_smoothed = smoothn(value_i, R2line+0.01);
%     
%     tx = rhos;
%     ty = value_smoothed;
%     
%     [tx, ia] = unique(tx);
%     ty = ty(ia);
%     s = spline(tx, ty, rho(bdi));    
%     bd_pos(i,1) =s;
%     
% end
% 
% % bd_pos= prune_bd_vis(bd_pos);
% 
% % smooth along the curve 
% 
% % % 
% % % x = bd_pos(:,1);
% % % xe = smoothn([x;x;x], [R2(bd);R2(bd);R2(bd)]);
% % % bd_pos(:,1) = xe(length(x)+1:2*length(x));
% % % 
% % % x = bd_pos(:,2);
% % % xe = smoothn([x;x;x], [R2(bd);R2(bd);R2(bd)]);
% % % bd_pos(:,2) = xe(length(x)+1:2*length(x));
% 
% 
% for i=1:size(bd_pos,1)
%     if norm(bd_vis0(i,:) - bd_pos(i,:))> changetol
%         bd_pos(i,:) = bd_vis0(i,:);
%     end
% end
% 
% 
% 
% end
% Smooth on the boundary  by line fit
function bd_pos = smooth_bd_vis_byfit( uv_pol, vis,  bd, R2)

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

% % prune outbound
% % prune bd_vis value to make sure it is in bd_vis_bd
% function bd_pruned = prune_bd_vis(bd_vis)
% bd_pruned = bd_vis;
% bd_pruned(bd_pruned(:,1)<0,1)=0;
% bd_pruned(bd_pruned(:,1)>8,1)=8;
% end

function mu_new = mu_chop(mu,bound,constant)
mu_new = mu;
ind = abs(mu)>=bound;
mu_new(ind) = constant*(cos(angle(mu(ind)))+1i*sin(angle(mu(ind))));
end


 