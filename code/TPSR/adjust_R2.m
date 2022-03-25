function R2_new = adjust_R2(face, vlabel, R2)
% reduce the vertex weight if it is near the boundary, because it may be
% shifted

vr = compute_vertex_ring(face);
R2_new = R2;
for i=1:length(vr)
   if sum(diff(vlabel(vr{i}) ))>1  % not all the same 
       
       R2_new(i) = R2(i)*0.1;
   else
       R2_new(i) = R2(i)*1.2;
   end
    
end 


end