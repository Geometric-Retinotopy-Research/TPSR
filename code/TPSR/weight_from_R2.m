%% weight_from_R2 function

function y = weight_from_R2(R) 
para_R0 =5; 
para_k = 0.1;
L =1;
y = L./(1+exp(-para_k*(R-para_R0))); 
y(R<para_R0)=0;
y(y<0)=0;
end