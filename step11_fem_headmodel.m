%%
clear;close all; clc;

%ID = {'004'; '005'; '007'; '008'; '009'; '012'; '013'; '014'; '015'; '016'; '017'; '019'; '020'; '021'; '022';}; %'001';'002'
ID = {'004'};

pathIn = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
pathOut = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
% pathIn = '/Volumes/NO NAME/';
% pathOut = '/Volumes/NO NAME/';
%%
for a = 1:size(ID,1)
    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_fem.mat']);
    
    cfg        = [];
    cfg.shift  = 0.3;
    cfg.method = 'hexahedral';
    mesh_fem = ft_prepare_mesh(cfg,mri_seg_fem);
    
    save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_fem_mesh'],'mesh_fem');
%%
    cfg = [];
    cfg.method ='simbio';
    cfg.conductivity = [1.79 0.33 0.43 0.01 0.14];   % order follows mesh.tissuelabel: csf gray scalp skull white
%%    
    headmodel_fem = ft_prepare_headmodel(cfg,mesh_fem);
    
    save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_fem_mesh_headmodel'],'headmodel_fem');

end