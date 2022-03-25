%%  This script load average result to create an atlas
clear;
addpath('freesurferReader');
addpath('utilities')
addpath(genpath('geometry-processing-package'));

%% load the osf pRF results
prfResults = load('../data/osfstorage-archive/prfresults.mat');
atlasResult = load('../data/osfstorage-archive/atlas.mat'); % manual drawn atlas
foveaid = 23683;
R2threshold = 5;
%% Show the first subject, of all data result.
lra = {'lh','rh'};
modelid =1;  % 1. Fit all of the available data (300 time points x 6 runs)  % 2. Fit the first half of each run (150 time points x 6 runs)% 3. Fit the second half of each run (150 time points x 6 runs)
for sid = 1:184   % sbuject 184 is all data used model
    
    fprintf('The %d-th subject name is %d\n',sid, prfResults.subjectids(sid));

    pRFSingleResult = prfResults.allresults(:,:,sid,modelid); 
    
    for looplr = 1:2
        
        lr = lra{looplr}; 
        % first read map info
        graycood2cortex=load('graycood2cortex');
        info = graycood2cortex.info;
        % read fs32k subject anatomical
        sub_folder = sprintf('../data/%d/MNINonLinear/fsaverage_LR32k',prfResults.subjectids(sid));
        
        mesh_sub = read_subject_32k_fs_LR(sub_folder);
        hemi_sub = getfield(mesh_sub, lr);
        
        
        
        
        
        
        
        hemi_sub.pRF = zeros(size(hemi_sub.inflated.vertices,1) ,6);
        hemi_sub.atlas_wang = zeros(size(hemi_sub.inflated.vertices,1) ,1);
        hemi_sub.atlas_hcp = zeros(size(hemi_sub.inflated.vertices,1) ,1);
        if(lr =='lh')
            hemi_sub.pRF(info.vertexInd{1}, 1:6) =   pRFSingleResult(1:info.vertexNum(1),:);
            hemi_sub.atlas_wang(info.vertexInd{1}) =  atlasResult.wang2015(1:info.vertexNum(1));
            hemi_sub.atlas_hcp(info.vertexInd{1}) =  atlasResult.glasser2016(1:info.vertexNum(1));
            
        else
            hemi_sub.pRF(info.vertexInd{2}, 1:6) =   pRFSingleResult(info.vertexNum(1)+1:info.vertexNum(1)+info.vertexNum(2),:);
            hemi_sub.atlas_wang(info.vertexInd{2}) =  atlasResult.wang2015(info.vertexNum(1)+1:info.vertexNum(1)+info.vertexNum(2));
            hemi_sub.atlas_hcp(info.vertexInd{2}) =  atlasResult.glasser2016(info.vertexNum(1)+1:info.vertexNum(1)+info.vertexNum(2));
            
        end
        
        
        color_ang_sub = prf_value_2_color(['polar_' lr], hemi_sub.pRF(:,1));
        color_ecc_sub = prf_value_2_color('ecc', hemi_sub.pRF(:,2));
        
        
        
        
        % write out 32k mesh
        atlascolor =atlasResult.glasser2016cmap(hemi_sub.atlas_hcp+1,:);
        
        % cut according to the inflated mesh
        
        hemi_roi = gf_set_roi_vertices_geodesic(hemi_sub, foveaid);
        
        
        F = double(hemi_roi.ROI.faces);
        V = double(hemi_roi.ROI.verticesPial);
        
        uv_cfm = disk_harmonic_map(F,V);
        
        uv_area = disk_omt(F,V,uv_cfm);
        
        Nv = size(hemi_roi.pial.vertices,1);
        uv = repmat([1, 1], [Nv, 1]);
        
        uv(hemi_roi.ROI.father,:) = uv_area;
        
        %
        
        
%         write_mfile([num2str(prfResults.subjectids(sid)) lr '_ecc.m'],'Face', hemi_sub.sphere.faces,...
%             'Vertex %d %f %f %f {uv=(%f %f) rgb=(%f %f %f) Vinflate=(%f %f %f) Vpial=(%f %f %f) Vsphere=(%f %f %f) prf=(%f %f %f %f %f %f) atlaswang=(%d) atlashcp=(%d)}\n', ...
%             [double(hemi_sub.inflated.vertices), uv , atlascolor,  double(hemi_sub.inflated.vertices), double(hemi_sub.pial.vertices), double(hemi_sub.sphere.vertices), hemi_sub.pRF(:,1:6) hemi_sub.atlas_wang hemi_sub.atlas_hcp]);
%         
%         write_mfile([num2str(prfResults.subjectids(sid))  lr '_ang.m'],'Face', hemi_sub.sphere.faces,...
%             'Vertex %d %f %f %f {uv=(%f %f) rgb=(%f %f %f) Vinflate=(%f %f %f) Vpial=(%f %f %f) Vsphere=(%f %f %f) prf=(%f %f %f %f %f %f)  atlaswang=(%d) atlashcp=(%d)}\n', ...
%             [double(hemi_sub.pial.vertices), uv, color_ang_sub,  double(hemi_sub.inflated.vertices), double(hemi_sub.pial.vertices), double(hemi_sub.sphere.vertices),hemi_sub.pRF(:,1:6)  hemi_sub.atlas_wang hemi_sub.atlas_hcp]);
%         
        
%         plot_surf(hemi_sub.inflated.faces, hemi_sub.inflated.vertices, color_ang_sub)
        
        
        
        %% lets show the atlas
        % uv_area = adj_area_map(F, V, uv_area0, [14545;16409], [0 0;0.5 0]);
        
        labels = hemi_roi.atlas_hcp(hemi_roi.ROI.father);
        
        interested = {'V1' 'V2' 'V3'  'V4'   };
        interested_id = [2 5 6 7 ]-1;   % Corrosponidng
        
        
        
        
        
        figure;
        R2= hemi_sub.pRF(hemi_roi.ROI.father, 5);
        plot_surf( hemi_roi.ROI.faces, uv_area, R2) ;
        color = color_ecc_sub(hemi_roi.ROI.father,:);
        color(R2<R2threshold,:) =nan;
        plot_surf( hemi_roi.ROI.faces, uv_area, color )
        hold on;
        Inid = [];
        for i = 1 : length(interested)
            
            loop = compute_region_boundary(F, labels, interested_id(i));
            
            Inid = [Inid; find(interested_id(i) ==labels)];
            plot_path(F,uv_area,loop,'g'); 
            
        end
        
        
        % cut the mesh
        if(sid>181)
        id2del  = [setdiff(1:length(labels), Inid) find(R2<R2threshold)'];
        else
            id2del  = [ find(R2<R2threshold)'];
        end
        [Fout, ~, vertex_father] = gf_remove_mesh_vertices( hemi_roi.ROI.faces,hemi_roi.ROI.verticesPial, id2del);
        
        atlas.faces = Fout;
        atlas.father = hemi_roi.ROI.father(vertex_father); % father relative to orginal mesh
        atlas.uv = uv_area(vertex_father,:);
        atlas.pRF = hemi_sub.pRF(atlas.father,:);
        atlas.verticesPial =   double(hemi_sub.pial.vertices(atlas.father,:)); %
        atlas.verticesInflated = double(hemi_sub.inflated.vertices(atlas.father,:)); %;
        atlas.verticesSphere =   double(hemi_sub.sphere.vertices(atlas.father,:));
        atlas.ang_color = prf_value_2_color(['polar_' lr], atlas.pRF (:,1));
        atlas.ecc_color = prf_value_2_color('ecc', atlas.pRF (:,2));
        write_mfile([num2str(prfResults.subjectids(sid)) lr '.m'],'Face', atlas.faces,...
            'Vertex %d %f %f %f {uv=(%f %f) rgb=(%f %f %f) Vinflate=(%f %f %f) Vpial=(%f %f %f) Vsphere=(%f %f %f) prf=(%f %f %f %f %f %f) atlaswang=(%d) atlashcp=(%d)}\n', ...
            [double(atlas.verticesInflated), atlas.uv, ...
             atlas.ang_color,...
            double(atlas.verticesInflated), double(atlas.verticesPial),...
            double(atlas.verticesSphere), atlas.pRF(:,1:6) hemi_sub.atlas_wang(atlas.father) hemi_sub.atlas_hcp(atlas.father)]);
        
        plot_surf(atlas.faces, atlas.uv, atlas.ang_color); title('Cut only V1 V2 V3 V4')
       
    end
    
end


%% Draw the ecc to radio
plot_surf(atlas.faces, atlas.uv, atlas.ang_color); title('Cut only V1 V2 V3 V4')
 
ang = atlas.pRF(:,1);
ecc = atlas.pRF(:,2);
vis_x = ecc.*cos(ang/180*pi);
vis_y = ecc.*sin(ang/180*pi);
% 
% landmarks =[];
% landpos =[] ;
%  for i=1:6
%     rid = find(abs(ecc - i)<0.5);
%     landmarks =[landmarks; rid];
%     dis = vecnorm(atlas.uv(rid,:)')';
%     landpos =[landpos ; atlas.uv(rid,1)./dis*i/8  atlas.uv(rid,2)./dis*i/8 ];
%  end
%  
%  uv_adj = disk_area_mapping(atlas.faces, atlas.verticesPial, atlas.uv, landmarks,landpos);
%  
figure
plot_surf(atlas.faces, atlas.uv, atlas.ecc_color); title(' |\mu| on V1 V2 V3 V4');
colorbar 
 