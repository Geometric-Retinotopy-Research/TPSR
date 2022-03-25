%%
% Description  -- function newid2keep = gf_expand_region(F, id2keep)
%       
%
% Parameter(s):
%		F[double array]  --  face of the region
%       id2keep[double array] -- current region
% 
% % return:
%       newid2keep[double array] -- new region after expand
%
%%
function newid2keep = gf_expand_region(F, id2keep)

vrs = compute_vertex_ring(F);

vre = [];
for i = 1: length(id2keep)
    vre = [vre vrs{id2keep(i)}];    
end

newid2keep = unique(vre);

end