%%
% Description  -- function [v1, v2]=logz(u,v, psnr_level)
%       log z model with gaussian noise
% Parameter(s):
%		u[double array] --  u in coordinate
%       v[double array] -- 
%       psnr_level[int] --
% return:
%       v1[double]  -- 
%       v2[double]  --
%% log z model with gaussian noise
function [v1, v2]=logz(u,v, psnr_level)

f= 0.5*log(u(:)+1i*v(:));    % add noise
v1 = real(f);
v2 = imag(f);

Rn = randn(length(v1),2);


if nargin < 3
    
else
    % add noise
    L = 1;
    while(1)
        v1n = v1 + L*Rn(:,1);
        v2n = v2 + L*Rn(:,2);
        if psnr(vecnorm([v1 v2]),vecnorm([v1n v2n])) < psnr_level
            L=L*0.9;
        else
            L=L*1.1;
        end
        
        v1n = v1 + L*Rn(:,1);
        v2n = v2 + L*Rn(:,2);
        if ( abs(psnr(vecnorm([v1 v2]),vecnorm([v1n v2n])) - psnr_level)  <0.01)
            L
            break;
        end
    end
    
    v1 = v1n;
    v2 = v2n;
     
end
 

if nargout < 2
    v1 =[v1 v2];
end
end
