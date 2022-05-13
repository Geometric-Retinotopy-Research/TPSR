%% prepare data
fn = '102311lh.m';
[Fm,Vm, Em]=read_mfile(['../data/mesh_data/' fn]);


[Ffull,Vfull, Efull]=read_mfile(['../data/mesh_data/' fn(1:end-2) '_ecc.m']);

uvm = disk_conformal_map(Fm,Vm);

roipatch = load('../data/v1');
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

%% topology smoothing
R2 = prf(:,5);
bd = compute_bd(Froi);
% anchorpos(:,1)=anchorpos(:,1);
% anchorpos(:,2)=anchorpos(:,2);
% visxy_s = topological_smooth_2(Froi,uv_p,visxy_corrected, R2, anchor, anchorpos);

changetol = 0.1; 
smooth_lambda0 = 0.001; % initial of smooth
smooth_avg_k = 2;  % size of smooth along the boundary
meanddth = 1; % average ddth 
visxy_s = topological_smoothing(Froi,uv_p,  visxy_corrected, R2,...
                                anchor, anchorpos,changetol, ...
                            smooth_lambda0,smooth_avg_k, meanddth);
						
%% plot the figure

figure
h=plot_surf(Froi,visxy_s,visxy_s(:,1)); hold on;
alpha(0.3)
h=plot_mesh(Froi,visxy_s,'Edgecolor','k'); hold on;
set(h,'FaceColor','none')

plot(visxy_s(anchor,1),visxy_s(anchor,2),'r-','Linewidth',2);
set(gca,'Fontsize',20);
set(gca, 'Color', 'none')

drawnow;
xlim([0 8])
ylim([1 5]);
