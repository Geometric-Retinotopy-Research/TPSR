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
%		uvw[double array]     --  uv coordinary
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
function [NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new, extra]= compute_fMRI_metrics(subjectfn, lr, hat_vis, uvw)

%% load pRF  pRF and fmri data
data_all = load_data(subjectfn, lr);

%% x y in pixel
[xhat, yhat]=vis2xy_pixel(hat_vis);


%% computer the map index


map2fmriid = uvw2dataid(uvw, data_all);
fn = fieldnames(data_all);
for k=1:numel(fn)
    if( isnumeric(data_all.(fn{k})) )
        if size(data_all.(fn{k}),1)>1             
            data.(fn{k}) = data_all.(fn{k})(map2fmriid,:);
        else
            data.(fn{k}) = data_all.(fn{k})(map2fmriid);
        end
    end
end
clear data_all



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




R2123 = [];
for vx = 1 : size(hat_vis,1)
   
    % For each run, collect the data and the model fit.  We project out polynomials
    % from both the data and the model fit.  This deals with the problem of
    % slow trends in the data.
    datats = {};
    modelts = {};
    for p=1:length(degs)
        datats{p} =  polymatrix{p}*data.fmri(vx,(p-1)*300+1:(p)*300)';
        parameters = [xhat(vx) yhat(vx) data.rfsize(vx)  data.gain(vx) data.expt(vx)];       
        modelts{p} = polymatrix{p}*modelfun(parameters,stimulusPP{p});
    end
    
    % Visualize the results  
    
     R2123(vx) = calccod(cat(1,modelts{:}), cat(1,datats{:}));     
     rmse123(vx) = rmse( (cat(1,modelts{:})+mdata(vx,:))./mdata(vx,:)', (cat(1,datats{:})+mdata(vx,:))./mdata(vx,:)');
     err123(vx,:) = cat(1,modelts{:}) - cat(1,datats{:});
      

    
end
 

% Check new R2
R2hat123 = [];
for vx = 1 : size(hat_vis,1)
    
    % For each run, collect the data and the model fit.  We project out polynomials
    % from both the data and the model fit.  This deals with the problem of
    % slow trends in the data.
    datats = {};
    modelts = {};
    
    rf_a = linspace(data.rfsize(vx)*0.1, data.rfsize(vx)*10, 50);

    for sss =1:length(rf_a)
        for p=1:length(degs)
        datats{p} =  polymatrix{p}*data.fmri(vx,(p-1)*300+1:(p)*300)';
        parameters = [xhat(vx) yhat(vx) rf_a(sss)  data.gain(vx) data.expt(vx)];       
        modelts{p} = polymatrix{p}*modelfun(parameters,stimulusPP{p});
        end
        tmp_R2 = calccod(cat(1,modelts{:}), cat(1,datats{:}));
        vgsss(sss) = tmp_R2;
    end
    
    % Visualize the results
   
      [R2hat123(vx),max_id] = max(vgsss);
      
      for p=1:length(degs)
          datats{p} =  polymatrix{p}*data.fmri(vx,(p-1)*300+1:(p)*300)';
          parameters = [xhat(vx) yhat(vx) rf_a(max_id)  data.gain(vx) data.expt(vx)];
          modelts{p} = polymatrix{p}*modelfun(parameters,stimulusPP{p});
      end
      rmsehat123(vx) = rmse( (cat(1,modelts{:}) + mdata(vx,:)')./ mdata(vx,:)', (cat(1,datats{:})+mdata(vx,:)')./mdata(vx,:)');       
      errhat123(vx,:) = cat(1,modelts{:}) - cat(1,datats{:});
end

 AIC_raw = calAIC(err123);
 AIC_new = calAIC(errhat123);
 
 R2_new=nanmean(R2hat123);
 R2_raw=nanmean(R2123);

 NRMSE_raw = nanmean(rmse123);
 NRMSE_new =  nanmean(rmsehat123);
  

extra.R2123 = R2123;
extra.R2hat123 = R2hat123;

extra.rmse123 = rmse123;
extra.rmsehat123 = rmsehat123;


extra.err123 = err123;
extra.errhat123 = errhat123;

end



% load 59k mesh and fMRI
function data = load_data(subjectfn, lr)
 
pRFSingleResult = load(['../data/prf_results_59k/' subjectfn '.mat']); % notice the results are in the unit of pixel

sid = str2double(subjectfn);


% read fsk subject anatomical 59k, because fmri is decoded in 59k mesh
% sub_folder = sprintf('../data/%d/MNINonLinear/fsaverage_LR59k',sid);



% push the subject to data  
% info = read_nii_info(sprintf('../data/%d/MNINonLinear/fsaverage_LR59k/%d.ArealDistortion_1.6mm_MSMAll.59k_fs_LR.dscalar.nii',sid,sid));
load('../data/graycood2cortex59k')
% mesh_sub = read_subject_59k_fs_LR(sub_folder); 
mesh_sub = pRFSingleResult.subject; 
if(lr =='lh')
    data.ang(info.vertexInd{1}) =  pRFSingleResult.results.ang(1:info.vertexNum(1));
    data.ecc(info.vertexInd{1}) =  pRFSingleResult.results.ecc(1:info.vertexNum(1));
    data.expt(info.vertexInd{1}) =  pRFSingleResult.results.expt(1:info.vertexNum(1));
    data.rfsize(info.vertexInd{1}) =  pRFSingleResult.results.rfsize(1:info.vertexNum(1));
    data.R2(info.vertexInd{1}) =  pRFSingleResult.results.R2(1:info.vertexNum(1));
    data.gain(info.vertexInd{1}) =  pRFSingleResult.results.gain(1:info.vertexNum(1));
    fmri =[];
    for i=1:6
        fmri = [fmri  double(pRFSingleResult.data{i}(1:info.vertexNum(1),:))];
    end
    data.fmri(info.vertexInd{1},:) = fmri;
    data.hemi  = mesh_sub.lh;
else
    data.ang(info.vertexInd{2}) =  pRFSingleResult.results.ang(info.vertexNum(1)+1:info.vertexNum(1)+info.vertexNum(2));
    data.ecc(info.vertexInd{2}) =  pRFSingleResult.results.ecc(info.vertexNum(1)+1:info.vertexNum(1)+info.vertexNum(2));
    data.expt(info.vertexInd{2}) =  pRFSingleResult.results.expt(info.vertexNum(1)+1:info.vertexNum(1)+info.vertexNum(2));
    data.rfsize(info.vertexInd{2}) =  pRFSingleResult.results.rfsize(info.vertexNum(1)+1:info.vertexNum(1)+info.vertexNum(2));
    data.R2(info.vertexInd{2}) =  pRFSingleResult.results.R2(info.vertexNum(1)+1:info.vertexNum(1)+info.vertexNum(2));
    data.gain(info.vertexInd{2}) =  pRFSingleResult.results.gain(info.vertexNum(1)+1:info.vertexNum(1)+info.vertexNum(2));
     fmri =[];
    for i=1:6
        fmri = [fmri  double(pRFSingleResult.data{i}(info.vertexNum(1)+1:info.vertexNum(1)+info.vertexNum(2),:))];
    end
    data.fmri(info.vertexInd{2},:) = fmri;
    data.hemi  = mesh_sub.rh;
end

data.sigma =  data.rfsize.*sqrt(data.expt);
end


%define a function map the mesh vertex id to fmri index
% use  individual subject space first read map info
function id = uvw2dataid(uvw32k, data)
Nq = size(uvw32k,1);
id = zeros(Nq,1);
for i=1:Nq
    dv = (data.hemi.sphere.vertices - uvw32k(i,:));
    [~,id(i)]=min(dot(dv,dv,2));
end
end



% recovery from hat_vis to xy in pixel
function [x,y] = vis2xy_pixel(hat_vis)

% from hat_vis to new x, y in visual image pixel
hat_ecc = hat_vis(:,1)/8*100;
hat_ang = hat_vis(:,2);

% recovery from extended
hat_ang_reverse = hat_ang;
hat_ang_reverse( hat_ang>pi/2 & hat_ang<3*pi/2) = hat_ang( hat_ang>pi/2 & hat_ang<3*pi/2)-pi; % V1
hat_ang_reverse( hat_ang>0 & hat_ang<pi/2) = -hat_ang( hat_ang>0 & hat_ang<pi/2); % V2d
hat_ang_reverse( hat_ang>-pi/2 & hat_ang<0) = hat_ang( hat_ang>-pi/2 & hat_ang<0); %V3d
hat_ang_reverse( hat_ang>3*pi/2 & hat_ang<2*pi) = -hat_ang( hat_ang>3*pi/2 & hat_ang<2*pi)+2*pi; %V2v
hat_ang_reverse( hat_ang>2*pi) = hat_ang( hat_ang>2*pi) - 2*pi; %V3v
hat_ang_reverse(hat_ang_reverse<0) = hat_ang_reverse(hat_ang_reverse<0) + 2*pi;

%
eccpixel =  hat_ecc;
anginrad =  hat_ang_reverse/180*pi;
    
x = eccpixel.*cos(anginrad) + 100.5;
y = 100.5 - eccpixel.*sin(anginrad);

end
