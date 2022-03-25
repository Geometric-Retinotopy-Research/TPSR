%% Set ROI in the struct of hemi
%% Syntax
  
%% Description

% hemi:  
 

%% Contribution
%  Author : Yanshuai Tu
%  Created: 2019/2/25
%  Revised: Not yet
%
%  Copyright @ Geometry Systems Laboratory
%  School of Computing, Informatics, and Decision Systems Engineering, ASU
%  http://gsl.lab.asu.edu/

%% Reference :
 
function ind2del = gf_get_del_id_geodesic(hemi,foveaid, radius)

if nargin <3
    radius =100;
end
    
start_points = foveaid;
vertex = hemi.pial.vertices;
faces = hemi.pial.faces;
[D,S,Q] = perform_fast_marching_mesh(double(vertex), double(faces), start_points);
 
% Loop thru hemispheres
ind2del  = find(D > radius);