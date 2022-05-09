clear;clc;close all;

addpath(genpath('geometry-processing-package'));
addpath('utilities');
addpath('topological-smoothing');

subjects = dir('../data/mesh_data/*lh.m');
rng(0)
%%
subi = 1;

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

%% Before Smoothing
[meanse, std_vd, mean_ang, std_ang, flip]=evaulate_metric(visxy_corrected,visxy_corrected,Froi,uv_p);
fprintf('Before smoothing meanse = %f, std_vd = %f, meanang = %f, maxang = %f flip =%d\n', meanse, std_vd, mean_ang, std_ang, flip);
 


vis_c = restore_vis(visxy_corrected);
showfigure(Froi, vis_c, anchor)
   
showresulton_inflatedmesh(Froi,Em.Vertex_Vinflate(vfather,:), uv_roi, vis_c, Ffull, Efull);


%% Average Smoother

vr = compute_vertex_ring(Froi);
visxy_s = [];
for i = 1:length(vr)
    visxy_s(i,:) = mean(visxy_corrected(vr{i},:));
end

[meanse, std_vd, mean_ang, std_ang, flip]=evaulate_metric(visxy_s,visxy_corrected,Froi,uv_p);
fprintf('Average smoothing meanse = %f, std_vd = %f, meanang = %f, maxang = %f flip =%d\n', meanse, std_vd, mean_ang, std_ang, flip);
 
vis_c = restore_vis(visxy_s);
showfigure(Froi, vis_c, anchor)   
showresulton_inflatedmesh(Froi,Em.Vertex_Vinflate(vfather,:), uv_roi, vis_c, Ffull, Efull);

%% Median

visxy_s = [];
for i = 1:length(vr)
    visxy_s(i,:) = median(visxy_corrected(vr{i},:));
end

[meanse, std_vd, mean_ang, std_ang, flip]=evaulate_metric(visxy_s,visxy_corrected,Froi,uv_p);
fprintf('Median smoothing meanse = %f, std_vd = %f, meanang = %f, maxang = %f flip =%d\n', meanse, std_vd, mean_ang, std_ang, flip);


vis_c = restore_vis(visxy_s);
showfigure(Froi, vis_c, anchor)   
showresulton_inflatedmesh(Froi,Em.Vertex_Vinflate(vfather,:), uv_roi, vis_c, Ffull, Efull);

%% Laplacian

visxy_s = laplacian_smooth(uv_p,visxy_corrected, prf(:,5));

[meanse, std_vd, mean_ang, std_ang, flip]=evaulate_metric(visxy_s,visxy_corrected,Froi,uv_p);
fprintf('Laplacian smoothing meanse = %f, std_vd = %f, meanang = %f, maxang = %f flip =%d\n', meanse, std_vd, mean_ang, std_ang, flip);
 
vis_c = restore_vis(visxy_s);
showfigure(Froi, vis_c, anchor)   
showresulton_inflatedmesh(Froi,Em.Vertex_Vinflate(vfather,:), uv_roi, vis_c, Ffull, Efull);
%%  Proposed

R2 = prf(:,5);
anchorid = compute_bd(Froi);


changetol = 0.5;
smooth_lambda0 = 0.002;
smooth_avg_k = 2;
meanddth = 1;
visxy_s = topological_smoothing(Froi,uv_p,  visxy_corrected, R2,...
    anchorid, anchorpos,changetol, ...
    smooth_lambda0,smooth_avg_k, meanddth);


[meanse, std_vd, mean_ang, std_ang, flip]=evaulate_metric(visxy_s,visxy_corrected,Froi,uv_p);
fprintf('Proposed smoothing meanse = %f, std_vd = %f, meanang = %f, maxang = %f flip =%d\n', meanse, std_vd, mean_ang, std_ang, flip);
 

vis_c = restore_vis(visxy_s);
showfigure(Froi, vis_c, anchor)   
showresulton_inflatedmesh(Froi,Em.Vertex_Vinflate(vfather,:), uv_roi, vis_c, Ffull, Efull);
 


 
showresulton_inflatedmesh(Froi,Em.Vertex_Vinflate(vfather,:), uv_roi, vis_c, Ffull, Efull);
view(90,0)
zoom(0.2) 

pos=get(gca,'Position');
pos(1) = pos(1) -0.2;
pos(2) = pos(2) +0.2;
set(gca, 'Position', pos)


 
%% show the inflated surface together with smoothed visual coordinate.
function showresulton_inflatedmesh(face, vinflated, uv, visxy, Ffull,Efull)
figure(123)
subplot(121)
r=[0.5:0.5:6];
theta=linspace(-pi/2-pi,pi/2+pi,30);
[c1,h1]=tricontour(face,uv(:,1),uv(:,2),visxy(:,1),r);
subplot(122)
[c2,h2]=tricontour(face,uv(:,1),uv(:,2),visxy(:,2),theta);


% construct a map from uv domain back to inflated mesh
for i=1:3
    Fitfun{i}=scatteredInterpolant(uv(:,1), uv(:,2), vinflated(:,i));
end

figure('Position', [100 100 400 500])
%  plot_surf(Ffull,Efull.Vertex_Vinflate,Efull.Vertex_rgb, 'Edgecolor','none', ...
%      'FaceColor',[189 140 124 ]/255,'FaceLighting','gouraud',...
%       'AmbientStrength',0.5);
   
rendercolor = Efull.Vertex_rgb;
rendercolor= rendercolor*0.8 + 0.2;
id = Efull.Vertex_atlaswang > 6 | Efull.Vertex_atlaswang<1;
repn = length( find(id));
rendercolor(id,:)= repmat([189 140 124 ]/255,repn,1);




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

view(16,-20)
axis off
zoom(3)




pos=get(gca,'Position');
pos(1) = pos(1) +0.2;
pos(2) = pos(2) -0.2;
set(gca, 'Position', pos)

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




function showfigure(F,visxy,anchor)
% compute the xy visual coordinate
  
visxy = [visxy(:,1)  visxy(:,2) ];

figure
h=plot_surf(F,visxy,visxy(:,1)); hold on;
alpha(0.2)
h=plot_mesh(F,visxy,'Edgecolor','k'); hold on;
set(h,'FaceColor','none')

plot(visxy(anchor,1),visxy(anchor,2),'r-','Linewidth',2);
set(gca,'Fontsize',20);
set(gca, 'Color', 'none')

drawnow; 
xlim([0,8])
ylim([-6,6])
end


 

function vis_c = restore_vis(visxy)
vis_c = [visxy(:,1) visxy(:,2)-pi];
end