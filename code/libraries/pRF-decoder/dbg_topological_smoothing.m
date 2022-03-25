 clear;clc;close all;

addpath(genpath('geometry-processing-package'));
addpath('utilities');
addpath('topological-smoothing');

subjects = dir('../data/mesh_data/*lh.m');
rng(0)
%%
subi = 2;

global closecontour
closecontour =1;
%% Prepare data
fn = subjects(subi).name;
[Fm,Vm, Em]=read_mfile(['../data/mesh_data/' fn]);

% Load the full hemisphere
[Ffull,Vfull, Efull]=read_mfile(['../data/mesh_data/' fn(1:end-2) '_ecc.m']);
uvm = disk_harmonic_map(Fm,Vm);

% load cut info

roipatch = load('../data/v1v2v3');
id2delete = roipatch.id2delete;


[Froi, V_roi, vfather] = gf_remove_mesh_vertices(Fm, Vm, id2delete);
uv_roi = uvm(vfather,:);
prf = Em.Vertex_prf(vfather,:);

visxy_corrected =correct_vis(Em, 'lh');
visxy_corrected = visxy_corrected(vfather,:);

[uv_p1, uv_p2] = cart2pol(uv_roi(:,1), uv_roi(:,2));
uv_p = [uv_p2, -uv_p1]; % r, theta
% anchor = roipatch.anchor;
anchor = compute_bd(Froi);
anchorpos = roipatch.anchorpos;


%% Porposed

R2 = prf(:,5);
anchorid = compute_bd(Froi);

% abnormalid = visxy_corrected(:,1) > 10; 
% visxy_corrected(abnormalid,:) = NaN;
% R2(abnormalid) = 0;

changetol = 0.5;
smooth_lambda0 = 0.002;
smooth_avg_k = 2;
meanddth = 1;


visxy_s = topological_smoothing_test2(Froi,uv_p,  visxy_corrected, R2, anchorid, anchorpos,changetol, smooth_lambda0,smooth_avg_k, meanddth);



[meanse, maxse, mean_ang, max_ang, flip]=evaulate_metric(visxy_s,visxy_corrected,Froi,uv_p);
fprintf('Proposed smoothing meanse = %f, maxse = %f, meanang = %f, maxang = %f flip =%d/%d\n', meanse, maxse, mean_ang, max_ang, flip,size(Froi,1));
 
vis_c = restore_vis(visxy_s);
showresulton_inflatedmesh(Froi,Em.Vertex_Vinflate(vfather,:), uv_roi, vis_c, Ffull, Efull);
camva(7.977609);
set(gca,'CameraPosition',[427.768368 -226.449049 -154.299598])
set(gca,'Position',[0.130000 0.110000 0.775000 0.815000])
set(gca,'Xlim',[-40.916702 -9.536835])
set(gca,'Ylim',[-103.025135 -24.170776])
set(gca,'Zlim',[-20.212266 35.609159])
set(gcf,'position',[680.000000 558.000000 560.000000 420.000000]) 
     
 
%% show the inflated surface together with smoothed visual coordinate.
function showresulton_inflatedmesh(face, vinflated, uv, visxy, Ffull,Efull)
figure(123)
subplot(121)
r=[0.5:0.5:6];
theta=-pi/2-pi:pi/12:pi/2+pi;
[c1,h1]=tricontour(face,uv(:,1),uv(:,2),visxy(:,1),r);
subplot(122)
[c2,h2]=tricontour(face,uv(:,1),uv(:,2),visxy(:,2),theta);


% construct a map from uv domain back to inflated mesh
for i=1:3
    Fitfun{i}=scatteredInterpolant(uv(:,1), uv(:,2), vinflated(:,i));
end

figure 
%  plot_surf(Ffull,Efull.Vertex_Vinflate,Efull.Vertex_rgb, 'Edgecolor','none', ...
%      'FaceColor',[189 140 124 ]/255,'FaceLighting','gouraud',...
%       'AmbientStrength',0.5);
   
rendercolor = Efull.Vertex_rgb*0;
repn = length( find(sum(rendercolor,2)==0));
rendercolor(sum(rendercolor,2)==0,:)= repmat([189 140 124 ]/255,repn,1);

% rendercolor =   [189 140 124 ]/255 ;
  plot_surf(Ffull,Efull.Vertex_Vinflate,'FaceVertexCData',rendercolor, 'Edgecolor','none','FaceLighting','gouraud',...
      'AmbientStrength',0.5);
  alpha(0.9)
%   
% light('Position',[500 0 0],'Style','local')
light('Position',[-500 0 0],'Style','local')
light('Position',[0 500 0],'Style','local')
light('Position',[0 -500 0],'Style','local')
light('Position',[0  0 -500 ],'Style','local')
light('Position',[0  0 500 ],'Style','local')
 
axis off 

 

axis manual
 hold on;
for i = 1:length(h1)
   Hi = h1(i) ;   
   plot3( Fitfun{1}(Hi.Vertices(:,1), Hi.Vertices(:,2)), ...
         Fitfun{2}(Hi.Vertices(:,1), Hi.Vertices(:,2)), ...
         Fitfun{3}(Hi.Vertices(:,1), Hi.Vertices(:,2)),'-','Linewidth',2,'color', [ 30  220 30]/255);
end

for i = 1:length(h2)
   Hi = h2(i) ;   
   plot3( Fitfun{1}(Hi.Vertices(:,1), Hi.Vertices(:,2)), ...
         Fitfun{2}(Hi.Vertices(:,1), Hi.Vertices(:,2)), ...
         Fitfun{3}(Hi.Vertices(:,1), Hi.Vertices(:,2)),'-','Linewidth',2,'color',[114 30 219]/255);
end

close(123); 

end

 
%% restore flip
function vis_c = restore_vis(visxy)
vis_c = [visxy(:,1) visxy(:,2)-pi];
end