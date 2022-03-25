%%
% Description  -- function [v_bd,k] = get_boundary_vertex(Fs,Es,i)%  
%       get the verticels of boundary
%
% Parameter(s):
%		Fs[double array]  --  
%		Es[double array]  --  
%		i[int]            --  
% 
% % return:
%       v_bd[double] -- verticles of boundary
%       k[int]  --  key index of boundary
%
%%
function [v_bd,k] = get_boundary_vertex(Fs,Es,i)%  
face_bnd_id = Es.Vertex_atlaswang(Fs);
    
face_bnd_type = sum(face_bnd_id==i,2); % 3: all inside  2 boundary 1 boundary  0 outside
bd_fid = find(face_bnd_type>1);
v_bd = 1/3*(Es.Vertex_uv(Fs(bd_fid,1),:) + Es.Vertex_uv(Fs(bd_fid,2),:) + Es.Vertex_uv(Fs(bd_fid,3),:));
k = boundary(v_bd(:,1),v_bd(:,2),1);
v_bd = v_bd(k,:);

end
