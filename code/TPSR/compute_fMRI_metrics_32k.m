%%
% Description  -- function [NRMSE_raw, NRMSE_new, ...
%           R2_raw_mean,R2_new_mean, ...
%           AIC_raw, AIC_new, ...
%           pc_raw, pc_new] = compute_fMRI_metrics(fMRI, pRF, pRF0)
%		Compute metrics by fMRI.
%
% Parameter(s):
%		subjectfn[String]     --  file name, read data from this file.
%		lr[String]            --  type of hemisphere.
%		hat_vis[double array] --  x,y in pixel.
%		roiid[int]            --  roi id.
% return:
%		NRMSE_raw[double]   --  Normalized RMSE of the pRF solution in fMRI	fitting.
%		NRMSE_new[double]   --  Normalized RMSE of the pRF0 solution in fMRI fitting
%		R2_raw[double] --  Average Valiance explained for pRF solution
%		R2_new[double] --  Average Valiance explained for pRF0 solution
%		AIC_raw[double]     --  AIC metric for pRF solution
%		AIC_new[double]     --  AIC metric for pRF0 solution
%		extra[double]       --  extra information.
%
%%
function [NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new, extra]= compute_fMRI_metrics_32k(subjectfn, lr, hat_vis, roiid)

%% load pRF  pRF and fmri data
data_all = load_data(subjectfn, lr);
fn = fieldnames(data_all);
for k=1:numel(fn)
    if( isnumeric(data_all.(fn{k})) )
        if size(data_all.(fn{k}),1)>1
            data.(fn{k}) = data_all.(fn{k})(roiid,:);
        else
            data.(fn{k}) = data_all.(fn{k})(roiid);
        end
    end
end
clear data_all

%% x y in pixel
[xhat, yhat]=vis2xy_pixel([hat_vis(:,1) hat_vis(:,2)-pi]); % the reason +pi is we shifted from it
[x, y]=vis2xy_pixel([data.ecc data.ang/180*pi]);




%% refit fMRI

% Define the model function.  This function takes parameters and stimuli as input and
% returns a predicted time-series as output.  Specifically, the variable <pp> is a vector
% of parameter values (1 x 5) and the variable <dd> is a matrix with the stimuli (frames x pixels).
% Although it looks complex, what the function does is pretty straightforward: construct a
% 2D Gaussian, crop it to <res>, compute the dot-product between the stimuli and the
% Gaussian, raise the result to an exponent, and then convolve the result with the HRF,
% taking care to not bleed over run boundaries.

stim_run_fns ={'RETCCW', 'RETCW','RETEXP','RETCON', 'RETBAR1', 'RETBAR2'};
for i=1:6
    stim_fn = sprintf('pRF-decoder/apertures/%ssmall.mat',stim_run_fns{i});
    stimstr= load(stim_fn);
    stimulus{i} = stimstr.stim;
end

TR =1; hrf = getcanonicalhrf(TR,TR)';          % HRF that was used in the model
degs =[3 3 3 3 3 3];  % vector of maximum polynomial degrees used in the model
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
    polymatrix{p} = projectionmatrix(constructpolynomialmatrix(300,0:degs(p)));
end

modelfun = @(pp,dd) conv2run(posrect(pp(4)) * (dd*[vflatten(placematrix(zeros(res),makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*abs(pp(3))^2))); 0]) .^ posrect(pp(5)),hrf,dd(:,prod(res)+1));



% Check old R2
for i=1:6
    mdata(:,(i-1)*300+1:(i)*300) = repmat(mean(data.fmri(:,(i-1)*300+1:(i)*300),2), [1,300]); % mean value respectively
end

    persistent parameters_array
    if isempty(parameters_array)
        for vx=1:size(hat_vis,1)
            parameters_array{vx} = NaN;
        end
    end    
  
    function [R2 , modelfit, fmridata, ...
            R2_new , modelfit_new]= R2by(x,y, xhat, yhat)        
        
        
        for i=1:6
            data_prf{i} = data.fmri(vx,(i-1)*300+1:(i)*300);
        end
        
        if isnan(parameters_array{vx})
            tr = 1;
            option.seedmode = [1];
            option.display = 'off';
            result = analyzePRF(stimulus,data_prf,tr, option);
            parai = squeeze(result.params);
            parameters_array{vx} = parai;           
        else
            parai = parameters_array{vx};           
        end
           
            
            datats = {};
            modelts = {};
            for pp=1:length(degs)
                parameters = [x, y, parai(3:5)]; % in pixel /8*100
                datats{pp} =  polymatrix{pp}*data.fmri(vx,(pp-1)*300+1:(pp)*300)';
                modelts{pp} = polymatrix{pp}*modelfun(parameters,stimulusPP{pp});
                
            end
        
        % Visualize the results
        modelfit = cat(1,modelts{:});
        modelfit = modelfit /  std(modelfit);
        fmridata = cat(1,datats{:});
        fmridata = fmridata / std(fmridata);
        R2 = calccod(modelfit, fmridata);
        if R2>100
            R2 =0;
        end
        
        
        datats = {};
        modelts = {};
        for pp=1:length(degs)
            parameters = [xhat, yhat, parai(3:5)]; % in pixel /8*100
            datats{pp} =  polymatrix{pp}*data.fmri(vx,(pp-1)*300+1:(pp)*300)';           
            modelts{pp} = polymatrix{pp}*modelfun(parameters,stimulusPP{pp});
           
        end
        % Visualize the results
        modelfit_new = cat(1,modelts{:});
        modelfit_new = modelfit_new /  std(modelfit_new); 
        R2_new = calccod(modelfit_new, fmridata);
         if R2_new>100
            R2_new =0;
        end
    end


R2data = [];
R2data_new =[];
rid = 1 : 10:size(hat_vis,1);
for vx = rid
    
    % For each run, collect the data and the model fit.  We project out polynomials
    % from both the data and the model fit.  This deals with the problem of
    % slow trends in the data.
    
    [R2 , modelopt, fmridataopt,...
        R2_new , modelfit_new] = R2by(x(vx), y(vx), xhat(vx),yhat(vx));
    R2data(vx) = R2; % negative is because we search minof -R
    rmsedata(vx) = rmse( (modelopt+mdata(vx,:))./mdata(vx,:)', (fmridataopt+mdata(vx,:))./mdata(vx,:)');
    errdata(vx,:) = modelopt - fmridataopt;
    R2data_new(vx) = R2_new; % negative is because we search minof -R
    rmsedata_new(vx) = rmse( (modelfit_new+mdata(vx,:))./mdata(vx,:)', (fmridataopt+mdata(vx,:))./mdata(vx,:)');
    errdata_new(vx,:) = modelfit_new - fmridataopt;
    
end



AIC_raw = calAIC(errdata(rid, :));
AIC_new = calAIC(errdata_new(rid, :));

R2_new=nanmean(R2data_new(rid));
R2_raw=nanmean(R2data(rid));

NRMSE_raw = nanmean(rmsedata(rid));
NRMSE_new =  nanmean(rmsedata_new(rid));


extra.R2123 = R2data;
extra.R2hat123 = R2data_new;

extra.rmse123 = rmsedata;
extra.rmsehat123 = rmsedata_new;


extra.err123 = errdata;
extra.errhat123 = errdata_new;

end


 
% load 32k mesh and fMRI
function data = load_data(subjectfn, lr)

pRFSingleResult = load(sprintf('F:/HCP32Kmat/%s_32k_pRF_fMRI.mat',subjectfn));
mesh_sub = pRFSingleResult.mesh_sub;

hemi = mesh_sub.(lr);
n = size(hemi.pRF,1);
data.ang = hemi.pRF(:,1);
data.ecc = hemi.pRF(:,2);
data.rfsize =  hemi.pRF(:,6);
data.expt = 0.5 + zeros(n,1);
data.gain = hemi.pRF(:,3);
data.R2 = hemi.pRF(:,5);
data.sigma =  data.rfsize.*sqrt(data.expt);
fmri = [];
for i=1:6
    fmri = [fmri  double(hemi.fMRI{i})];
end

data.fmri = fmri;
end



% recovery from hat_vis to xy in pixel
function [x,y] = vis2xy_pixel(hat_vis)

% from hat_vis to new x, y in visual image pixel
hat_ecc = hat_vis(:,1)/8*100;
hat_ang = hat_vis(:,2);
 
eccpixel =  hat_ecc; 
y = eccpixel.*cos(hat_ang) + 100.5;
x = 100.5 - eccpixel.*sin(hat_ang);

end
