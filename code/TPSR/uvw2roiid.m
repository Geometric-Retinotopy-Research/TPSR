%%
function roiid = uvw2roiid(V,uvw)
roiid = zeros(size(uvw,1),1);
for i=1:size(uvw,1)
    [mindis, id]= min(vecnorm(V' - uvw(i,:)'));
    roiid(i) = id;
end
end
