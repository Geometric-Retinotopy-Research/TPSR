%%
% Description  -- function loop = compute_region_boundary(F,labels, r_label)
%					compute a loop for the given mesh, with region formed by  r_label == labels
%					The idea is to cut out r_label ~= labels vertices, then compute_boundary%
% Parameter(s):
%		F[double array]       --  array of faces.
%		labels[double array]  --  labels of the faces.
%		r_label[double array] --  r_labels of the faces.
% return:
%       loop[double array]    --  return a loop for the given mesh, with region formed by  r_label == labels
%%

function loop = compute_region_boundary(F,labels, r_label)




Flabel = labels(F);

fid2keep = sum( double(Flabel == r_label),2)==3;

Fregion = F(fid2keep,:);

loop = compute_bd(double(Fregion));
 
end


 
%%
%  
% function  test_code
% %%
% F = hemi_sub.inflated.faces;
% V = hemi_sub.inflated.vertices;
% labels = hemi_sub.atlas_hcp;
% 
% for i = 1 : 181
%     if atlasResult.glasser2016labels{i}(1)=='V' || atlasResult.glasser2016labels{i}(1)=='v'
%         atlasResult.glasser2016labels{i}
%       i
%         loop = compute_region_boundary(F, labels, i-1);
%         plot_path(F,V,loop,'r');
%         hold on;
% break;
%     end
%     
% end
%  %%
% 
% end