%%
% Description  -- function Vi = fsaverage2fslr(Vq, lr)
%       convert fsaverage sphere coordinate to fs_LR
%
% Parameter(s):
%		Vq[String]            --  Vq
%		lr[String]            --  type of hemisphere.
% 
% return:
%       Vi[double] -- return sphere coordinate for lr
%
%%
function Vi = fsaverage2fslr(Vq, lr)
if lr == 'lh'
    giifile=gifti('standard_mesh_atlas/L.sphere.32k_fs_LR.surf.gii');
    Vlr = double(giifile.vertices);
    
    giifile=gifti('standard_mesh_atlas/resample_fsaverage/fs_LR-deformed_to-fsaverage.L.sphere.32k_fs_LR.surf.gii');
    Vavg = double(giifile.vertices);
    
else
    if   lr == 'rh'
        giifile=gifti('standard_mesh_atlas/R.sphere.32k_fs_LR.surf.gii');
        Vlr = double(giifile.vertices);
        
        giifile=gifti('standard_mesh_atlas/resample_fsaverage/fs_LR-deformed_to-fsaverage.R.sphere.32k_fs_LR.surf.gii');
        Vavg = double(giifile.vertices);
    else
        error('only lr or rh allowed');
        
    end
end

Vi = sphere_interpretation(Vavg, Vlr, Vq);