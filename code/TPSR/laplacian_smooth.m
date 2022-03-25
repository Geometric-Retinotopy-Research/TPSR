%%
% Description  -- function [v_bd,k] = get_boundary_vertex(Fs,Es,i)%  
%       get the verticels of boundary
%
% Parameter(s):
%		uv[double array]  --  uv coordinate
%		f[double array]   --  face of the region
%		W[int]            --  optional parameter
% 
% % return:
%       v_bd[double] -- verticles of boundary
%       k[int]  --  key index of boundary
%
%%
function fn = laplacian_smooth(uv,f, W)
if ~exist('W','var')
    W = ones(size(f,1),1);
end

Mres = 200; Nres = 200;


x = {linspace(min(uv(:,1)),max(uv(:,1)),Mres),linspace(min(uv(:,2)),max(uv(:,2)),Nres)};
[xx,yy] = ndgrid(x{1},x{2});
grid_vec = [xx(:) yy(:)];
id0 = xx.^2+yy.^2 >1; % take r>1 points



Fw = scatteredInterpolant(uv(:,1), uv(:,2), W);
W = reshape(Fw(grid_vec),[Mres,Nres]);
W(id0)=0;
W(W<0)=0;
    
for i = 1:size(f,2)
    % Create   data on grid, ecc,R2, ang
    feature = f(:,i);
    Fr = scatteredInterpolant(uv(:,1), uv(:,2), feature);
    f_on_grid = reshape(Fr(grid_vec),[Mres,Nres]);
%     f_on_grid(id0)=0;
    
    %  Do smoothing
    smooth_f_on_grid = smoothn(f_on_grid,W);
    fn(:,i) = griddata(xx, yy, smooth_f_on_grid, uv(:,1), uv(:,2));
end

end