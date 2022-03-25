%% Set ROI in the struct of hemi
%% Syntax
 
% hemi = gf_set_roi_vertices(hemi,Y)
%% Description

% hemi:  hemi sphere, contains inflated and sphere and pial; cut the mesh
% by geodesic 
 

%% Contribution
%  Author : Yanshuai Tu
%  Created: 2019/2/25
%  Revised: Not yet
%
%  Copyright @ Geometry Systems Laboratory
%  School of Computing, Informatics, and Decision Systems Engineering, ASU
%  http://gsl.lab.asu.edu/

%% Reference :
 
function roi_out = gf_remove_del_id(hemi_in,ind2del)
 

% First, remove the vertex from mesh.
F = hemi_in.inflated.faces;
V = hemi_in.pial.vertices;
Nv = size(V,1);
[Fout, ~, vertex_father] = gf_remove_mesh_vertices(F, V, ind2del);


ind2keep = setdiff( 1:Nv, ind2del);
verticesPial =     double(hemi_in.pial.vertices(ind2keep,:)); %

roi_out.F = Fout;
roi_out.V = verticesPial;

roi_out.father = vertex_father; %
roi_out.verticesPial = roi_out.V;

roi_out.verticesInflated =     double(hemi_in.inflated.vertices(ind2keep,:)); %

roi_out.verticesVery_Inflated =     double(hemi_in.very_inflated.vertices(ind2keep,:)); %

roi_out.verticesSphere =     double(hemi_in.sphere.vertices(ind2keep,:)); %  

roi_out.pRF = double(hemi_in.pRF(ind2keep,:));

roi_out.thickness = double(hemi_in.thickness(ind2keep,:));  

roi_out.curvature = double(hemi_in.curvature(ind2keep,:));

roi_out.aparc = double(hemi_in.aparc(ind2keep,:));

roi_out.atlas_wang = double(hemi_in.atlas_wang(ind2keep,:));

roi_out.atlas_hcp = double(hemi_in.atlas_hcp(ind2keep,:));
  
         

roi_out.description ='V is equivalent to  verticesPial for convience ';
