%% Description  -- function uv = disk_conformal_map(face, vertex, roi)
%     nake a conformal map from a disk plot
%
%% parameter(s): 
%      face[double array]  -- connectivity of mesh
%      vertex[double array]  -- vertex of mesh	
%% return: 
%      uv[double array]  -- coordinate of mesh
% 
%% 
%disk conformal map: if roi is given, the vertex in roi is refined
function uv = disk_conformal_map(face, vertex, roi)
uv = disk_harmonic_map(face,vertex);  
 
mu =  compute_bc(face, uv, vertex) ; 
bd = compute_bd(face);
bd_uv = uv(bd,:);

for i=1:20

    uv =linear_beltrami_solver(face, uv, mu, bd, bd_uv);
    mu = compute_bc(face, uv, vertex) ; 

end


if nargin == 3
    
    roi_f = find(ismember(face(:,1), roi) & ismember(face(:,2), roi) & ismember(face(:,3), roi));
%     fprintf('mean mu = %f\n', nanmean(abs(mu(roi_f))));
    for i=1:20
        
        mu = compute_bc(face, uv, vertex) ;        
        mun = 0*(mu);
        mun(roi_f) = mu(roi_f);           
        uv0 = uv;
        uv =linear_beltrami_solver(face, uv, mun, bd, bd_uv);        
        id = isnan(uv(:,1));
        if ~isempty(id)
            uv = uv0;
            break
        end
%         mu = compute_bc(face, uv, vertex) ;        
%         fprintf('mean mu = %f\n', mean(abs(mu(roi_f)))); 
        
    end
    
end



end
 
