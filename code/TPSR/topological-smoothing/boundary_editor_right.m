%% boundary editot:
% Generate the solve boundary, and a fix reasonable value, with ecc and
% polar value

close all;clc;clear
addpath(genpath('../'))


global para_R0 para_k
para_R0 =25;
para_k = 0.1;


fn = '999999rh.m';
[Fm,Vm, Em]=read_mfile(['../data/mesh_data/' fn]);

uvm = disk_conformal_map(Fm, Vm);
%
% uvm = disk_area_mapping(Fm, Vm, uv_c, [4154], [0,0]);



uv = Em.Vertex_uv;
prf = Em.Vertex_prf;

R2 = prf(:,5);


face = Fm;
vertex = uv;
%
% figure
% plot_surf(face,vertex,  prf_value_2_color('ecc', prf(:,2)) )
figure
% plot_surf(face,vertex,   )
plot_surf(face,vertex,  prf_value_2_color('rh', prf(:,1))*0.2+0.8* Em.Vertex_rgb)

axis off
xy =[];
while 1
    
    [x,y, button] = ginput(1);
    xy = [xy; x y];
    hold on ;
    
    plot(xy(:,1),xy(:,2),'w-')
    
    if(button==3)
        break
    end
    
end

% fn = input('filename:','s');
fn = 'v1v2v3right'
%% infer the boundary retinotopic  values

x = Em.Vertex_uv(:,1);
y = Em.Vertex_uv(:,2);
[in, on] = inpolygon(x,y,xy(:,1),xy(:,2));
Nv = size(Vm,1);
id2delete = setdiff(1:Nv, find(in));
[F, V, vfather] = gf_remove_mesh_vertices(Fm, Vm, id2delete);
bd = compute_bd(F);
[in, on] = inpolygon(x,y,x(vfather(bd)),y(vfather(bd)));
id2delete = setdiff(1:Nv, find(in));

uv = uvm(vfather,:);

prf = Em.Vertex_prf(vfather,:);



figure
plot_mesh(F, uv, 'FaceVertexCData', Em.Vertex_rgb(vfather,:))
colorbar
title('boundary ')


% infer the boundary value

%%
bd = compute_bd(F);

anchor = vfather(bd);

anchors =[anchor;anchor(1)];

figure
plot_surf(face,vertex, prf_value_2_color('rh',  Em.Vertex_prf(:,1)) )
hold on;
plot(Em.Vertex_uv(anchors,1), Em.Vertex_uv(anchors,2),'w-')


[theta,rho]=cart2pol(Em.Vertex_uv(:,1), Em.Vertex_uv(:,2));

Em.Vertex_prf_copy = Em.Vertex_prf;

%% show vis
vis= correct_vis(Em, 'rh') ;
plot_mesh(F, [uv, vis(vfather,2)  ], 'FaceVertexCData', Em.Vertex_rgb(vfather,:))
axis auto
% hold on; plot_mesh(F, [uv, vis(vfather,2).* (Em.Vertex_atlaswang(vfather,:)==3) ], 'FaceVertexCData', Em.Vertex_rgb(vfather,:))


%%

for i=1:length(anchor)
    
    [thetai,rhoi]=cart2pol(Em.Vertex_uv(anchor(i),1), Em.Vertex_uv(anchor(i),2));
    
    id = find(abs(theta - thetai)<0.05);
    
    
    [thetas, sid ]= sort(rho(id));
    eccs = vis(id(sid),1);
    
    
    
    ang_smoothed = smoothn(eccs, weight_from_R2(R2(id(sid))));
    
    figure(102); hold off;
    plot(thetas, eccs, thetas, ang_smoothed); hold on;
    s = spline(thetas,ang_smoothed,rhoi);
    plot(rhoi, s, 'go')
    ylim([0, 10])
    % fit a function
    
    anchorpos(i,1) =s;
    
end
%
for i=1:length(anchor)
    
    [thetai,rhoi]=cart2pol(Em.Vertex_uv(anchor(i),1), Em.Vertex_uv(anchor(i),2));
    
    id = find(abs(rho - rhoi)<0.05);
    
    
    [thetas, sid ]= sort(theta(id));
    angs = vis(id(sid),2);
    
    
    
    R2line = 20*weight_from_R2(R2(id(sid)));
    
    R2line(abs(diff(angs))>0.5)=0;
    
    mm = movmean(R2line,8);
    R2line(mm<0.9)=0;
    
    
    ang_smoothed = smoothn(angs, R2line);
    
    figure(102); hold off;
    plot(thetas, angs, thetas, ang_smoothed); hold on;
    s = spline(thetas,ang_smoothed,thetai);
    plot(thetai, s, 'go')
    plot(thetas,0.1*R2line,'b-');
    % fit a function
    
    anchorpos(i,2) =s;
    
end

%%
n =length(anchor);
anchorposext = [anchorpos;anchorpos;anchorpos];
anchorpossmooth= [smooth(anchorposext(:,1),3) smooth(anchorposext(:,2),3)];
anchorpos =  anchorpossmooth(n+1:2*n,:);
%


save(['../' fn], 'id2delete','anchor',   'anchorpos')


figure
%
plot(anchorpos(:,1), anchorpos(:,2))
%%

