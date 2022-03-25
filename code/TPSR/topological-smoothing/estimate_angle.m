%%
% Description  -- function [hat_vis, flip, smooth_lambda] = topological_smoothing(face, uv_pol, vis, R2, bd_id, bd_vis0, changetol, smooth_lambda0, smooth_avg_k, meanddth)
%       Estimate the angle of two tangent vectors

%
% Parameter(s):
%		F[double array] -- faces of himesphere
%       V[double array] -- vertex of himesphere
%       vis[double array] -- Correct the visual coordinates
%       gridx[double array] -- x of grid
%       gridy[double array] -- y of grid
% 
% return:
%       angs[double array] -- angle after smoothing processing
%		vispos[int] -- visual coordinates
%
%%
% Estimate the angle of two tangent vectors
function [angs,vispos] = estimate_angle(F,V, vis, gridx, gridy)
% estimate angle by 2 vertex ring

figure (200); clf
%  plot_surf(F, V, prf_value_2_color('ecc',vis(:,1)))); hold on; 
alpha(0.5)
[c1,~]=tricontour(F,V(:,1),V(:,2),vis(:,1),gridx);
axis equal; 
i = 1 ;
j=1;
H1 = (struct);
while 1
    id = i+1:i+c1(2, i);
    H1(j).Vertices = c1(:,id)';
    H1(j).Level = c1(1,i);
    
    i = i+1+c1(2,i);
    j = j + 1;
    
    if i>= size(c1,2)
        break
    end
end

C1 = (struct);
for i=1:length(H1)
    
    xx = H1(i).Vertices(:,1);
    nanid = isnan(xx);
    xxfit=xx(~nanid);
    yy = H1(i).Vertices(:,2);
    yyfit=yy(~nanid);
    
    try
        
        tfit = (1:length(xxfit))';
        fx=fit(tfit,xxfit,'smoothingspline');
    
        fy=fit(tfit,yyfit,'smoothingspline'); 
        
        C1(i).fx = fx;
        C1(i).fy = fy;
        C1(i).level =  H1(i).Level;
        C1(i).span = 0:0.5:length(xxfit);       
        
    catch
        C1(i).level  = NaN; 
    end
end


[c2,~]=tricontour(F,V(:,1),V(:,2),vis(:,2),gridy);
axis equal; hold on;


i = 1 ;
j=1;
H2 = (struct);
while 1
    id = i+1:i+c2(2, i);
    H2(j).Vertices = c2(:,id)';
    H2(j).Level = c2(1,i);
    
    i = i+1+c2(2,i);
    j = j + 1;
    
    if i>= size(c2,2)
        break
    end
end

C2 = (struct);
for i=1:length(H2)
    
    
    xx = H2(i).Vertices(:,1);
    nanid = isnan(xx);
    xxfit=xx(~nanid);
    yy = H2(i).Vertices(:,2);
    yyfit=yy(~nanid);
    
    try
        tfit = (1:length(xxfit))';
        fx=fit(tfit,xxfit,'smoothingspline');
    
        fy=fit(tfit,yyfit,'smoothingspline'); 
        C2(i).fx = fx;
        C2(i).fy = fy;
        C2(i).level =  H2(i).Level;
        C2(i).span = 0:0.5:length(xxfit);
        
    catch
        C2(i).level  = NaN; 
    end
end

angs = zeros( length(C1)*length(C2),1);
vispos = zeros( length(C1)*length(C2),2);
i=1;
for c1 = C1
    for c2 = C2
        if isnan(c1.level) || isnan(c2.level)
            angs(i) =NaN;
        else
            
        [pos,v1,v2]=get_intersection_pos_tangent(c1,c2); 
        angs(i) = acos(dot(v1,v2)/sqrt(dot(v1,v1)*dot(v2,v2)));
        vispos(i,:) = [c1.level, c2.level];
        end
        i = i+1;
        
    end
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

end