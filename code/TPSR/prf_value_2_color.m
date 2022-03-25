function color = prf_value_2_color(vartype, value, maxvalue)
%%
% color = prf_value_2_color(vartype, value)
%
% Calculate the pseudo-colour for value.
%   
%   vartype :     'Polar_Lh', 'Polar_Rh', 'Eccen', 'Sigma', or 'R^2'
%   value :       angles in rad!!
%
% Yanshuai:  Modified from samsrf_colourcode
% 
% Avoid outbound by setting value>maxvalue as maxvalue
 

% Avoid nan, we set nan to black
nanid = isnan(value);
value(isnan(value))=0;

if(~exist('maxvalue'))
    maxvalue = max(value(:));    
end


if strcmpi(vartype, 'polar_lh') || strcmpi(vartype, 'lh')
    t = value;    
    t = t + 270;  % calibrate angles to stimulus starting position
    t = mod(ceil(t),360) + 1;   % ensure between 1-360
    cmap = (fscol);
    color = cmap(t,:);
elseif strcmpi(vartype, 'polar_rh')|| strcmpi(vartype, 'rh')
    t = value;    
    t = -t - 270;  % calibrate angles to stimulus starting position
    t = mod(ceil(t),360) + 1;   % ensure between 1-360
    cmap = (fscol);
    color = cmap(t,:);
elseif strcmpi(vartype, 'ecc')
    
    value(value<0)=0;
    value(value>maxvalue)=maxvalue;
    nr = ceil(value/maxvalue*360); % normalized rho in degrees
    nr(nr==0) = 1;
    
    cmap = (fscol);
    cmap = [cmap(315:360,:); cmap(1:314,:)];
    color = cmap(nr,:);
else
    
   error('Not supported yet!\n') ;
end

color(nanid,:)= color(nanid,:)*0;
end
