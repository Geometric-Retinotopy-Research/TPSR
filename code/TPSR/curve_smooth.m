%%
% Description  -- function smoothed = curve_smooth(anchorpos, R2)
%		smoothed the curve
% Parameter(s):
%       anchorpos[double array]	  --  anchor position
%       lr[String]                --  type of the hemisphere.
% return:
%		smoothed[double array]    --  results of smoothed
%
%%
function smoothed = curve_smooth(anchorpos, R2)
R2(R2<10) =0;
R2 = R2+min(R2);

anchorposdense = dense_arc(anchorpos);
w = dense_arc(R2);
n =size(anchorposdense,1);
anchorposnew = anchorposdense;

adj =1;
for t =1:200
    
    [wn,seg] = winding_number(anchorposnew);
    if wn <1
        % out if good
        break
    end
    
    % fix winding number
    segd =[];
    for i=1:adj
        segd =[segd seg+i seg-i];
    end
    segd = unique(segd);    
    segd(segd>n) = n;
    segd(segd<1) = 1;
    
    w(segd) = w(segd)*0.1;    
    if (nanmean(w(segd))<1e-9)
        adj = adj+1;
        w = dense_arc(R2);
    end
       
    
    we = [w; w; w];
    we(we<0) =0;
    anchorpos_ext = [anchorposdense; anchorposdense;anchorposdense];
    anchorpos_ext_smooth = [smoothn(anchorpos_ext(:,1), we) smoothn(anchorpos_ext(:,2), we)];
    anchorposnew = anchorpos_ext_smooth(n+1:2*n,:);
            
end

smoothed = sparse_arc(anchorposnew);



figure
plot(anchorpos(:,1),anchorpos(:,2)); hold on;
plot(smoothed(:,1),smoothed(:,2));

end

function dpos = dense_arc(pos)

% dense point
n = size(pos,1);
dpos = zeros(2*n,size(pos,2));
for i=1:size(pos,2)
    dpos(:,i) = spline(1:2:2*n,pos(:,i), 1:2*n);
end

end


function spos = sparse_arc(pos)

% dense point
n = size(pos,1);
spos = zeros(n/2,size(pos,2));
for i=1:size(pos,2)
    spos(:,i) = spline(1:n,pos(:,i), 1:2:n);
end

end