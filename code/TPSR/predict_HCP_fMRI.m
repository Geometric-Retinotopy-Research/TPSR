%%
% Description  -- function yt = predict_HCP_fMRI(x,y,s, g, expt)
%       compute the fMRI based on x, y, and sigma x,y, in the unit of pixel
% 
% Parameter(s):
%		x[int]  --  
%       y[int]  --  
%       s[int]  -- sigma
%       g[int]  -- 
%       expt[int]  --
% return:
%       yt[double array]  -- result of predict
%
%%

function yt = predict_HCP_fMRI(x,y,s, g, expt)
 
% Calculate energy based frmi and para
persistent stimulus
if isempty(stimulus)
    stim_run_fns ={'RETCCW', 'RETCW','RETEXP','RETCON', 'RETBAR1', 'RETBAR2'};
    for i=1:6
        stim_fn = sprintf('pRF-decoder/apertures/%ssmall.mat',stim_run_fns{i});
        stimstr= load(stim_fn);
        stimulus{i} = stimstr.stim;
    end
end

% Prepare the stimuli for use in the model
persistent stimulusPP
if isempty(stimulusPP)
    stimulusPP = {};
    for p=1:length(stimulus)
        stimulusPP{p} = squish(stimulus{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
        stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
    end
end


persistent hrf
TR =1;
if isempty(hrf)
    hrf = getcanonicalhrf(TR,TR)';          % HRF that was used in the model
end


degs =[3 3 3 3 3 3];  % vector of maximum polynomial degrees used in the model
res = size(stimulus{1}(:,:,1));
resmx = max(res);

persistent xx yy
if isempty(xx) || isempty(yy)
    % Pre-compute cache for faster execution
    [~,xx,yy] = makegaussian2d(resmx,2,2,2,2);
end


% Construct projection matrices that fit and remove the polynomials.
% Note that a separate projection matrix is constructed for each run.
persistent polymatrix
if isempty(polymatrix)
    polymatrix = {};
    for p=1:length(degs)
        polymatrix{p} = projectionmatrix(constructpolynomialmatrix(300,0:degs(p)));
    end
end
persistent modelfun
if isempty(modelfun)
%     modelfun = @(pp,dd) conv2run((dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]),hrf,dd(:,prod(res)+1));
    modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));

end


 
% For each run, collect the data and the model fit.  We project out polynomials
% from both the data and the model fit.  This deals with the problem of
% slow trends in the data.

modelts = {};
pp = [x, y, s, g, expt];
for p=1:length(degs)
    modelts{p} = polymatrix{p}*modelfun(pp,stimulusPP{p});
end
yt = cat(1,modelts{:});

end


