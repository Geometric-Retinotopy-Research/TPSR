%%		
% Description  -- function nmse = calcNRMSE(y,x)
%		calculating NRMSE value.
%
% Parameter(s): 
%		y[double]  --  
%		x[double]  --  
% return: 
%		nmse[double]  -- return back the NRMSE value.
% 
%		
%%
function nmse = calcNRMSE(y,x)
	nmse=mse(x,y)/mean(abs(x));
end