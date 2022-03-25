% compute the flip triangles
function [flip, flipid] = get_flips(face, uv, fn)
flipid = find(abs(compute_bc(face,uv,fn))>1);
flip = length(flipid);
end