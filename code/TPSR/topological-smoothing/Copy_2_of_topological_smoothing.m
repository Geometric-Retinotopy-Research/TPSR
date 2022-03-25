
function [hat_vis, flip, smooth_lambda] = topological_smoothing(face, uv_pol, vis, R2, bd_id, bd_vis0, changetol, smooth_lambda0, smooth_avg_k, meanddth)
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
%% update boundary value by predict from interier

if changetol > 0 
    for ti = 1:3
        bd_vis_new = smooth_bd_vis_byfit( uv_pol, vis,  bd_id, R2);
        for i=1:size(bd_vis_new,1)
            bd_vid_change = bd_vis_new(i,:) - bd_vis0(i,:);
            if norm(bd_vid_change)> changetol
                bd_vis_new(i,:) = bd_vis0(i,:) + bd_vid_change/norm(bd_vid_change)*changetol; % update to the maxi tolerance
            end
        end
        bd_vis = bd_vis_new;
    end
else
    
    bd_vis = bd_vis0;
end
%% topology smoothing
smooth_lambda = smooth_lambda0; 
fixed = 0;
for iteration = 1:20

    % check topology condition
    flip = get_flips(face, uv_pol, hat_vis);
    if flip > 0  
        % smooth and chop if not topological
        [hat_vis, flip] = make_diffeomorphic(face, uv_pol, hat_vis, bd_id, bd_vis);
        
        if size(hat_vis,1) ~= size(vis,1)
            hat_vis = vis;
        end
    end
    
    % recheck the topology condition    
    if( flip > 0)  
        
        [~,flipid] = get_flips(face, uv_pol, hat_vis);        
        % fail to fix topology, the only possible reason is boundary 
        % value is not good, so we smooth on boundary value
                
        % smooth on boundary locally
        flipid_v = [face(flipid,1); face(flipid,2);face(flipid,3)];        
        % feed flip vertices, and smooth on boundary near flip
        bd_vis_new = smooth_bd_vis_bybd(bd_id, bd_vis,smooth_avg_k, flipid_v);
        for i=1:size(bd_vis,1)
            bd_vid_change = bd_vis_new(i,:) - bd_vis(i,:);
            if norm(bd_vid_change)> changetol
                bd_vis_new(i,:) = bd_vis(i,:) + bd_vid_change/norm(bd_vid_change)*changetol; % update to the maxi tolerance
            end
        end
        bd_vis = bd_vis_new;
%         figure(108); hold off;
%         plot_path(face,hat_vis, bd_id); hold on;
%         plot(bd_vis(:,1),bd_vis(:,2),'g-');
        
        % increase the strenth of smooth
        if  changetol >= 0
            smooth_avg_k = smooth_avg_k +1;
        end
    end
   
     
          
    % laplacian smooth with the given boundary condition
    A = (diag(R2) + smooth_lambda* laplace_beltrami(face,uv_pol));
    b = (diag(R2)*hat_vis);
    in = setdiff(1:size(uv_pol,1), bd_id);
    hat_vis(in,:) = A(:,in)\(b - A(:,bd_id)* bd_vis);
    hat_vis(bd_id,:) = bd_vis;
%     for i=1:size(uv_pol,1)
%         
%         if norm(hat_vis(i,:) - vis(i,:))>1  && R2(i)>30
%             hat_vis(i,:) = 0.2*(hat_vis(i,:) - vis(i,:))/norm(hat_vis(i,:) - vis(i,:)) + vis(i,:);
%         end
%     end   
 
    % Break the loop if we fixed topology and value is within tolerance
    [meandev, meandd] = visual_deviation(hat_vis, vis);      
    flip = get_flips_nobd(face, uv_pol, hat_vis);
    if( (flip < 1)  &&  meandd < meanddth)         
        fixed =1;
        smooth_lambda = smooth_lambda*1.01;        
        fprintf('Fixed topology, and the deviation is  %f\n', meandd);     
        break;
    end
    
    % dispy how much deviation from the raw data
    flip = get_flips_nobd(face, uv_pol, hat_vis);
    [meandev,mean_dd] = visual_deviation(hat_vis, vis);    
%     fprintf('Iteration %d: mean|smooth_f -f0| = %f #%d of overalaps \n',iteration,mean_dd, flip);
    
    
end

% do it again with bigger boundary deviation
if ~fixed
    try
        close(108)
    catch
    end

    changetol = changetol*1.5; % such tol is not big enough.  
    smooth_lambda0 = smooth_lambda0*1.5;    
    meanddth = meanddth*1.5; 
    if smooth_avg_k >= 120
        fprintf('Failed')
        return        
    else
        [hat_vis, flip, smooth_lambda] = topological_smoothing(face, uv_pol, vis, R2, bd_id, bd_vis, changetol, smooth_lambda0, smooth_avg_k, meanddth);        
    end
 
end

end


% make_diffeomorphic:  try to make diffeomorphic by smooth and chop mu
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
fv = face_area(face,fn); 
fu = face_area(face,uv); 
face(ismember(face(:,1),bd)| ismember(face(:,2),bd) |ismember(face(:,3),bd) | fv<1e-4| fu<1e-4,: )=[]; 
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



smoothsz = 5 + round(kspan/20);

for i = 1:n
    bdi = bd(i);
    if ismember(bdi, flipidv)
        % if flip related boundary smooth locally
        bd_visext(n + i-smoothsz:i+ smoothsz + n,1) = smooth(bd_visext(n + i-smoothsz:i+ smoothsz + n,1), kspan);
        bd_visext(n + i-smoothsz:i+ smoothsz + n,2) = smooth(bd_visext(n + i-smoothsz:i+ smoothsz + n,2), kspan);
        
    else
        continue;
    end
end

bd_pos = bd_visext(n + 1:2*n,:);

end
 

function mu_new = mu_chop(mu,bound,constant)
mu_new = mu;
ind = abs(mu)>=bound;
mu_new(ind) = constant*(cos(angle(mu(ind)))+1i*sin(angle(mu(ind))));
end


 