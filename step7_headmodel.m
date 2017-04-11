clear;close all; clc;

ID = {'002';} % '003'; '004'; '005'; '007'; '008'; '009'; '012'; '013'; '014'; '015'; '016'; '017'; '019'; '020'; '021'; '022';}; %'001';'002'


%pathIn = '/Volumes/OUTPUT DATA/Sian/EEG_coords_reconstruction/';
%pathOut = '/Volumes/OUTPUT DATA/Sian/EEG_coords_reconstruction/';
pathIn = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
pathOut = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';

for a = 1:size(ID,1)

    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_mesh']);

    %% Prepare head model (volume conduction model)
    cfg = [];
    cfg.method = 'bemcp';
    headmodel_bem = ft_prepare_headmodel(cfg,mesh_bem);
    
    save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_mesh_bem'],'headmodel_bem');
    
end