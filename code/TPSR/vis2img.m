%% Description  --  value2image(F,V, vis, rect, fn)
%		convert values to image, and save image to file path.		
% 
% Parameter(s): 
%     F[double array]     -- connectivity of mesh
%     V[double array]     -- vertex of mesh
%	  vis[int]            -- v value
%     rect[double array]  -- location and size infomation
%     fn[string]          -- file name to save
%% 
function img_box = vis2img(vis,R2)

lb = min(vis);
ru = max(vis);

M = 500;



% 
% [X,Y]=meshgrid(linspace(lb(1),ru(1),M), linspace(lb(2),ru(2),M));
% 
 R2(R2<0)=0;
% V=scatteredInterpolant(vis,R2/max(R2),'natural');
% 
imgbox.img = zeros(M,M,1);
imgbox.lb = lb;
imgbox.ru = ru;

for i=1:size(vis,1)
    rx = round((vis(i,1) - lb(1))/(ru(1)-lb(1))*M) ;
    ry = round((vis(i,2) - lb(2))/(ru(2)-lb(2))*M) ;
    if rx <1 
        rx =1;
    end
    if rx >M 
        rx =M;
    end
    
    if ry <1 
        ry =1;
    end
    if ry >M 
        ry =M;
    end
    
    imgbox.img(rx,ry) = R2(i);
end


figure;
imshow(imgbox.img')

end