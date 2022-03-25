%%
% Description  -- function [R2, yd, yf]= evaulate_analyzeprf_fitness(stimulus,data, tr, parameters)
%		compute a result of prf analyze
% Parameter(s):
%       stimulus[double array]	  --  stimulus to retina
%       data[douhle array]        --  
%		parameters[double array]  --  
% 
% return:
%		R2[double array]  --  
%       yd[double array]  --
%       yf[double array]  --
% 
%%
function [R2, yd, yf]= evaulate_analyzeprf_fitness(stimulus,data, tr, parameters)
if ~iscell(stimulus)
  stimulus = {stimulus};
end
hrf = getcanonicalhrf(tr,tr)';
for i=1:length(stimulus)
    degs(i) = round(size(stimulus{1},3)*tr/60/2);
end
res = size(stimulus{1}(:,:,1));
resmx = max(res);
% Pre-compute cache for faster execution
[~,xx,yy] = makegaussian2d(resmx,2,2,2,2);

% Prepare the stimuli for use in the model
stimulusPP = {};
for p=1:length(stimulus)
    stimulusPP{p} = squish(stimulus{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
    stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
end

% Construct projection matrices that fit and remove the polynomials.
% Note that a separate projection matrix is constructed for each run.
polymatrix = {};
for p=1:length(degs)
    Ts = size(stimulus{p},3);
    polymatrix{p} = projectionmatrix(constructpolynomialmatrix(Ts,0:degs(p)));
end

model = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));
 

R2 = [];
yf = [];
yd = [];
for vx = 1 : size(data,1)    
    % For each run, collect the data and the model fit.  We project out polynomials
    % from both the data and the model fit.  This deals with the problem of
    % slow trends in the data.
    datats = {};
    modelts = {}; 
    for p=1:length(degs)
        Ts = size(stimulus{p},3);
        datats{p} =  polymatrix{p}*data(vx,(p-1)*Ts+1:(p)*Ts)'; 
        modelts{p} = polymatrix{p}*model(parameters(:,vx),stimulusPP{p});
    end    
    
    datacat = cat(1,datats{:});      
    modelfit = cat(1,modelts{:}); 
    modelfit = modelfit /  std(modelfit);
    
   
    
    % Visualize the results    
    yf(:, vx) = modelfit;
    yd(:, vx) = datacat;  
    
    R2(vx,1) = calccod(modelfit, datacat);     
end 


