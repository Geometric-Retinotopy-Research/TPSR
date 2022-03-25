%%
% Description  -- function ft = resample_parametric(uv,f, uvt)
%		
% Parameter(s):
%		uv[double array]  --  uv coordinate value
%       f[double array]  --  face 
%       uvt[int]  -- uv transpose
% 
% return:
%       ft[double array]  -- interpolant of the result
%
%%
function ft = resample_parametric(uv,f, uvt)
% uv Nx3 on sphere
% f  Nx1 
% return Nx1 ft

% transpose if not Nx3 but 3XN
if(size(uv,1) ==3 && size(uv,1) ~= 3 )
    uv = uv';
end

if(size(uvt,1) ==3 && size(uvt,1) ~= 3 )
    uvt = uvt';
end

x = uv(:,1);
y = uv(:,2);
z = uv(:,3);


xt = uvt(:,1);
yt = uvt(:,2);
zt = uvt(:,3);

[az,el, ~] = cart2sph(x,y,z);
[azt,elt,~] = cart2sph(xt,yt,zt);


ft=interp_mesh(az,el,f,azt,elt);

end


%% test




function V=interp_mesh(x,y,v,qx,qy)
    % Construct the interpolant:
    Fi = scatteredInterpolant(x,y,v);
    V=Fi(qx,qy);
    
end
