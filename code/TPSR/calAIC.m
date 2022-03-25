%%		
% Description  -- function aic = calAIC(err)
%		calculating AIC with error.
%
% Parameter(s): 
%		err[double]  --  error value which use calculate AIC.
%
% return: 
%		aic[double]  -- return back the AIC value.
% 
%		
%%
function aic = calAIC(err)

err = double(err);

% SS = sum(err(:).*err(:));
% N = length(err(:));
% K = size(err,1)*4;
% 
% aic = AIC(SS,N,K);


f = 0;
for i=1:size(err,1)
    
    mu = nanmean(err(i,:));
    sigma = nanstd(err(i,:));       
    for j =1:size(err,2)
        logg = log(g(err(i,j), mu, sigma));
        if abs(logg)<1e9
            f = f + logg;
        end
       
    end
end
k = size(err,1)*4;
aic = -2*f + 2*k;

aic = aic/size(err,1)/6;% we have six runs in each fitting

end

function y=g(x, mu, sigma)
y =  1/sigma/sqrt(2*pi)*exp(-0.5*(x-mu).^2/sigma^2);
if(isnan(y))
    y = 0;
end
end


function aic = AIC(SS,N,K)

%AIC Akaike's Information Criterion for curve-fitting model comparison
%   AIC(SS,N,K) where SS is the sum of squared residuals (as returned by e.g. lsqcurvefit)
%   N is the number of data points and K is the number of coefficients. 
%   Computes and returns the corrected AIC score for the model.
%   
%   The model with the lowest AIC value is the best fit.
%
%   References: (1) Motulsky, H. and Christopoulos, A. (2003) Fitting models to biological data using 
%   linear and nonlinear regression. A practical guide to curve fitting. Graphpad Software Inc., San Diego CA.
%   www.graphpad.com
%
%   (2) Akaike, H. (1974) "A new look at the statistical model identification." IEEE Transactions on Automatic Control, AC-19, 716-723
%
%   (3) Hurvich, C. M., and Tsai, C-L. (1989). Regression and time series model selection in small samples. Biometrika, 76, 297-307.
%   [the AIC correction]   
%
%   NOTE: this computes AIC from sum-of-squares (SS), and thus uses SS as
%   an esimator for the maximum likelihood; when actually fitting
%   distributions using MLE, then use AICL instead!
%
%   Mark Humphries 11/10/2004

K  = K + 1; % additional degree-of-freedom is model!

% raw AIC
aic = N .* log(SS./N) + 2 .* K;

% apply correction in case N close to K
aic = aic + (2.*K.*(K+1)) ./ (N - K - 1);


end