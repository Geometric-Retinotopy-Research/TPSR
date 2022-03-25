%% Description  -- disk_conformal_map(face,vertex)
%     nake a conformal map from a disk plot
%
%% parameter(s): 
%      face[double array]  -- connectivity of mesh
%      vertex[double array]  -- vertex of mesh	
%% return: 
%      uv[double array]  -- coordinate of mesh
% 
%% 
% Estimate the angle of two tangent vectors
function ang = estimate_angle_by_pos(F,V, vis, pos)
% estimate angle by 2 vertex ring

figure (200)
for i = 1:size(pos,1)
    
    try
    [c1,~]=tricontour(F,V(:,1),V(:,2),vis(:,1),[pos(i,1) pos(i,1)]); 
    catch
         ang(i) =NaN;
        continue
    end
    
    if(isempty(c1))
        ang(i) =NaN;
        continue
    end
        
    xx = c1(1,2:end);
    yy = c1(2,2:end);
    nanid = isnan(xx);
    xxfit=xx(~nanid);
    yyfit=yy(~nanid);    
    try
        
        tfit = (1:length(xxfit))';
        fx=fit(tfit,xxfit','smoothingspline');    
        fy=fit(tfit,yyfit','smoothingspline');         
        cv1.span = tfit;
        cv1.fx =fx;
        cv1.fy=fy;
    catch
        ang(i) =NaN;
        continue
    end
    
    try
    [c2,~]=tricontour(F,V(:,1),V(:,2),vis(:,1),[pos(i,2) pos(i,2)]); 
    catch
         ang(i) =NaN;
        continue
    end
    if(isempty(c2))
        ang(i) =NaN;
        continue
    end
    xx = c2(1,2:end);
    yy = c2(2,2:end);
    nanid = isnan(xx);
    xxfit=xx(~nanid);
    yyfit=yy(~nanid);    
    try
        
        tfit = (1:length(xxfit))';
        fx=fit(tfit,xxfit','smoothingspline');    
        fy=fit(tfit,yyfit','smoothingspline');         
        cv2.span = tfit;
        cv2.fx =fx;
        cv2.fy=fy;
    catch
        ang(i) =NaN;
        continue
    end
    
    [~,v1,v2]=get_intersection_pos_tangent(cv1,cv2);
    ang(i) = acos(dot(v1,v2)/sqrt(dot(v1,v1)*dot(v2,v2)));
end
 
global closecontour
if closecontour
    close(200)
end

end



function [pos,v1,v2]=get_intersection_pos_tangent(c1,c2)
[c1grid,c2grid]=meshgrid(c1.span, c2.span);
dis = sqrt((c1.fx(c1grid)-  c2.fx(c2grid)).^2 + (c1.fy(c1grid)-  c2.fy(c2grid)).^2 );
dis = reshape(dis,size(c1grid));
[min_val,idx]=min(dis(:));
if (min_val <0.1)
    [row,col]=ind2sub(size(dis),idx);
    pos = [c1.fx(c1grid(row,col)) c1.fy(c1grid(row,col))];
    v1 = [differentiate(c1.fx, c1grid(row,col)) differentiate(c1.fy, c1grid(row,col))];
    v2 = [differentiate(c2.fx, c2grid(row,col)) differentiate(c2.fy, c2grid(row,col))];
else
    pos =[nan nan];
    v1 =[nan nan];
    v2 =[nan nan];
end
v1 = v1/norm(v1);
v2 = v2/norm(v2);

end