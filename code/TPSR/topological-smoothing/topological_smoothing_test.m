% topological_smoothing:
% Inputs: 
% face: Face list
% uv_pol: parametric coordinate in polar representation
% vis: noisy vis coordinates
% R2: variance from fMRI fitting
% bd_id: boundary points
% bd_vis0: the initial boundary visual coordinate
% bd_tol: boundary change tolerance
% lambda0 : inital smooth parameter
% bd_avg_k: average along boundary 
% avg_dev_thre : overall average smooth threshold 

% Outputs:  
% hat_vis: the smoothed visual coordinates
% flip: the flip number after smooth, if failed

function [hat_vis, flip] = topological_smoothing_test(face, uv_pol, vis, R2, bd_id, bd_vis0, ...
                                                      bd_tol, lambda0,smooth_avg_k, avg_dev_thre)

% default parameters
if ~exist('avg_dev_thre','var') 
    avg_dev_thre = 1.0; 
end
if ~exist('smooth_avg_k','var')        
    smooth_avg_k = 2;
end
if ~exist('lambda0','var')
    lambda0 = 0.01;
end
if ~exist('bd_tol','var')
    bd_tol = 0;
end
  
abnormalid = vis(:,1) > 10; 
vis(abnormalid,:) = NaN;
R2(abnormalid) = -100;

 
%  W = weight_from_R2(R2); 

for to = 1:10
      W = weight_from_R2(R2); 
       hat_vis= vis;
     % make boundary simple
     fprintf('Fixing Boundary\n')
     for tb=1:10
         [hat_vis_fix, fixed_bd] = make_simple_boundary(face, uv_pol,hat_vis, W, bd_id, bd_tol );
         if fixed_bd
             hat_vis = hat_vis_fix;
             break;
         else
             fprintf('We cannot make such bd_tol, increase by 0.5\n');
             bd_tol = bd_tol*1.5;
         end
     end
     
     % smooth and fix topology interior
     fprintf('Fixing interior\n')
     for ti=1:3
         [hat_vis, fixed, mu_v] = make_topological_interior(face, uv_pol,hat_vis, W, bd_id);
         if fixed
             break
         else
             fprintf('Lower the weight of non-topological interior\n');
             W(mu_v > 0.99) = W(mu_v > 0.99)* 0.8;  % degrade the flipping points
             W = laplacian_smooth(uv_pol,W);
         end
     end
     
     % Good, they are topological, furthur check or optimize     
     vis_dis = visual_deviation(hat_vis, vis);
     vis_dis
     if  vis_dis < avg_dev_thre
         fprintf('Fixed topology, vis distortion=%f (Tolerance %f)\n',...
             visual_deviation(hat_vis, vis), avg_dev_thre);
         break
         
%      else
%          % move boundary          
%          
%           bd_lin = smooth_bd_vis_byfit( uv_pol, vis, bd_id, R2);
%           
%           hat_vis(bd_id,:) = (hat_vis(bd_id,:) + bd_lin)/2;
%           W(bd_id) =1;
%          
     end
       
 end
 
 
 figure
 plot_mesh(face, hat_vis)
 hold on
 plot(hat_vis(bd_id,1), hat_vis(bd_id,2),'-r')
 
 
 
 
end 
%%
function [hat_vis_fix, fixed, mu_v] = make_topological_interior(face, uv_pol,...
                                        hat_vis, W, bd_id)
fixed = 0 ; 
Operator = createOperator(face,uv_pol);
 hat_vis0  = hat_vis;
 mu_v =[];

for i=1:5
    
  
   
    % laplacian smooth with the given boundary condition
    hat_vis = laplacian_smooth(uv_pol, hat_vis, W);    
    hat_vis(bd_id,:) = hat_vis0(bd_id,:);
    
    % fix topology
    mu = compute_bc(face,uv_pol, hat_vis);
    mu_new = mu_chop(mu, 0.99, 0.9); 
    
    [hat_vis_fix,mu_new] = linear_beltrami_solver(face, uv_pol, mu_new,bd_id, hat_vis(bd_id,:));    
    flip = get_flips_nobd(face, uv_pol, hat_vis);     
    fprintf('\tInteration %d with %d flips\n',i,flip)
    if flip<1
        fixed = 1;
        fprintf('Fixed topology\n')
        break
    end
    
    
    mu_v = Operator.f2v *mu_new;

end


end

%%
function [hat_vis, fixed_bd] = make_simple_boundary(face, uv_pol, ...
                                                    hat_vis, W, bd_id, bd_tol )
fixed_bd = 0;
Operator = createOperator(face,uv_pol);
Wb = W(bd_id);
for i=1:5
    hat_vis = laplacian_smooth(uv_pol,hat_vis, W);
    mu = compute_bc(face, uv_pol, hat_vis);
    vmu = Operator.f2v*mu;
    % process boundary: make it a simple looop
    Wb = Wb.*decay_mu(vmu(bd_id));
    bd_vis = boundary_smooth(hat_vis(bd_id,:),  Wb);
    bd_vis = rect_vis(bd_vis, hat_vis(bd_id,:), bd_tol);
    hat_vis(bd_id,:) =  bd_vis;
    
    [n,seg]=winding_number(bd_vis);
    if  n ~= 0 
        Wb(seg(:)) = Wb(seg(:)) * 0.95;        
    else
        fixed_bd =1;
        fprintf('Fixed Boundary\n');
        break
    end
    
   
end

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

%% boundary smooth

function bd_vis = boundary_smooth(bd_vis0, w_bd)
dv = bd_vis0(2:end,:) - bd_vis0(1:end-1,:);
dl = sqrt(dot(dv,dv,2));
s = [0; cumsum(dl)/sum(dl)];

dl_e1 = sqrt(dot(bd_vis0(end,:) - bd_vis0(1,:), bd_vis0(end,:) - bd_vis0(1,:), 2));
ss = [s;s+dl_e1 + 1; s+ 2*dl_e1 + 2]; 
bd_vis = bd_vis0;
n = length(s);
for i=1:size(bd_vis0,2)
    yi = [bd_vis0(:,i); bd_vis0(:,i);bd_vis0(:,i)];
    yi_s = smooth1d(ss, yi, w_bd);
    bd_vis(:,i) = yi_s(n+1:2*n); 
end 
end

%% decay the weight by \mu
function d = decay_mu(mu)
d = exp(-abs(mu).^2/0.5);
end


%% 1D Lap smooth
function ys = smooth1d(x,y, w_b)
xq = linspace(0,3,1e3)';
yq = spline(x,y,xq); 
wq = spline(x,[w_b;w_b;w_b],xq);
wq(wq<0)=0;
ysq = smoothn(yq,wq); 
ys = spline(xq,ysq,x); 
end


%% estimate the visual_deviation, in the unit of visual degree
function [meanvd, mean_dd] = visual_deviation(hat_vis, vis)
vis_distance = sqrt( (hat_vis(:,1).* (hat_vis(:,2) - vis(:,2) )).^2 + (hat_vis(:,1) - vis(:,1) ).^2);
meanvd = nanmean(vis_distance);
mean_dd = nanmean(vecnorm( (hat_vis-vis)'));
end

%% compute the number of flips without bd
function flip = get_flips_nobd(face, uv, fn)
bd = compute_bd(face);
face(ismember(face(:,1),bd)| ismember(face(:,2),bd) |ismember(face(:,3),bd),:)=[]; 
flip = length(find(abs(compute_bc(face,uv,fn))>1));
end


%% compute the number of flips
function [flip, flipid] = get_flips(face, uv, fn)
flipid = find(abs(compute_bc(face,uv,fn))>1);
flip = length(flipid);
end

%%  Smooth on the boundary locally
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



%% Smooth on the boundary  by line fit
function bd_pos = smooth_bd_vis_byfit( uv_pol, vis, bd, R2)

rho = uv_pol(:,1);
theta = uv_pol(:,2);
bd_pos = zeros(length(bd),2);

R2(isnan(R2))=-5;
% for angle
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
    
    
    id = find( isfinite(value_i) & isfinite(thetas) & isfinite(R2line));
    
    
    % smooth   
    cf = fit(thetas(id),value_i(id),fittype('poly3'),'Weight',R2line(id)+0.01); 
    s = cf(uv_pol(bdi,2));    
    bd_pos(i,2) =s;
end


% for eccentricity
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
    id = find( isfinite(value_i) & isfinite(rhos) & isfinite(R2line));
    
    cf = fit(rhos(id),value_i(id),fittype('poly3'),'Weight',R2line(id)+0.01); 
   
    bd_pos(i,1) = cf(uv_pol(bdi,1));    
        
end


end

%% mu_chop(mu,bound,constant)
function mu_new = mu_chop(mu,bound,constant)
mu_new = mu;
ind = abs(mu)>=bound;
mu_new(ind) = constant*(cos(angle(mu(ind)))+1i*sin(angle(mu(ind))));
end


 