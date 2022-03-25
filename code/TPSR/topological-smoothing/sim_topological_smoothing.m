
function [hat_vis, flip] = topological_smoothing(face, uv, vis, R2, bd_id, bd_vis,changetol, smooth_lambda0, smooth_avg_k, meanddth)
%% initilize hat_v 
hat_vis = vis;     
%% topology smoothing
W = R2;
smooth_lambda = 1;
in_id = setdiff(1:size(uv,1), bd_id);

    % laplacian_smoothing with pointwise R
    function hat_vis = laplacian_smoothing(hat_vis, W)
        A = (diag(W) + smooth_lambda* laplace_beltrami(face,uv));
        b = diag(W)*hat_vis;
        hat_vis(in_id,:) = A(:,in_id)\(b - A(:,bd_id)* bd_vis);
        hat_vis(bd_id,:) = bd_vis;    
    end

    % topology_projection:  try to make diffeomorphic by smooth and chop mu
    function [hat_vis, flip] = topology_projection(hat_vis)        
        flip = get_flips(face, uv, hat_vis);
        if flip<1
            return
        end        
        oldhat = hat_vis;       
        upperBound = 0.95;       
        chopVal = 0.95;  % Chop on mu
        for i=1:1000
            % compute mu
            mu = compute_bc(face,uv, hat_vis); mu(isnan(mu))=0;            
            choped_mu = mu_chop(mu,upperBound,chopVal); 
            [lbs_vis,~] = linear_beltrami_solver(face,uv,choped_mu, bd_id, bd_vis);
            hat_vis = lbs_vis;
            % reverse back if nan occurs
            id = [find(isnan(hat_vis(:,1))); find(isnan(hat_vis(:,2)))];
            hat_vis(id,:) = oldhat(id,:) ;            
            flip = get_flips(face, uv, hat_vis);
            if( flip < 1 )
                break  % break when diffeomorphic
            end
            chopVal = chopVal-0.0001;
        end        
    end

    
for iteration = 1:20
    
    
    % laplacian smooth with the given boundary condition
    vis_smooth = laplacian_smoothing(hat_vis, W);
    
    % topology_projection
    hat_vis = topology_projection(vis_smooth);
    
    dvis = vecnorm(hat_vis'-vis_smooth')';
    dvis = (dvis-min(dvis))/(max(dvis)-min(dvis));
    W = W./(dvis+0.01);
    W = W * mean(R2)/mean(W);   
    
    
    % Break the loop if we fixed topology and value is within tolerance       
    flip = get_flips(face, uv, hat_vis);
    [~, vis_dis] = visual_deviation(hat_vis, vis);
    fprintf('Fixed topology, the deviation is  %f\n', vis_dis);
end
 
 

end


function [meanvd, mean_dd] = visual_deviation(hat_vis, vis)
vis_distance =sqrt( (hat_vis(:,1).* (hat_vis(:,2) - vis(:,2) )).^2 + (hat_vis(:,1) - vis(:,1) ).^2);
meanvd = mean(vis_distance);
mean_dd = mean(vecnorm( (hat_vis-vis)'));
end

% compute the number of flips without consider the boundary 
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
 

function mu_new = mu_chop(mu,bound,constant)
mu_new = mu;
ind = abs(mu)>=bound;
mu_new(ind) = constant*(cos(angle(mu(ind)))+1i*sin(angle(mu(ind))));
end


 