function  boundary_optomization(F123, uv_p, vis123, R2123,anchor, anchorpos)
global F123  uv_p  vis123  R2123 anchor  
x0 =[anchorpos(1:10:end,1); anchorpos(1:10:end,2)];
options = optimoptions('fminunc');
fun = @fitfun;
[x, fval] = fminunc(fun,x0,options)


end


function fitval  = fitfun(x)
global F123  uv_p  vis123  R2123 anchor  
warning('off')

anchorpos(:,1) = spline(1:10:length(anchor), x(1:length(x)/2), 1:length(anchor));
anchorpos(:,2) = spline(1:10:length(anchor), x(length(x)/2+1:end), 1:length(anchor));


% anchorpos(:,1) = smooth(anchorpos(:,1));
% anchorpos(:,2) = smooth(anchorpos(:,2));
% anchorpos = [x(1:length(anchor)) x(length(anchor)+1:end) ];
[hatvis,flip] = topological_smooth(F123, uv_p, vis123, R2123,anchor,anchorpos);

df = mean(dot(hatvis - vis123, hatvis - vis123, 2));
fitval = df + flip;
warning('on')
end