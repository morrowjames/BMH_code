clear;close all; clc;

%ID = {'004'; '005'; '007'; '008'; '009'; '012'; '013'; '014'; '015'; '016'; '017'; '019'; '020'; '021'; '022';}; %'001';'002'
ID = {'012';}; 

pathIn = '/Volumes/OUTPUT DATA/Sian/EEG_coords_reconstruction/';
pathOut = '/Volumes/OUTPUT DATA/Sian/EEG_coords_reconstruction/';
% pathIn = '/Volumes/NO NAME/';
% pathOut = '/Volumes/NO NAME/';

for a = 1:size(ID,1)


    %% Segment the anatomical data for FEM headmodel
    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice.mat']);

    cfg           = [];
    cfg.output    = {'gray','white','csf','skull','scalp'};
    cfg.brainsmooth = 2;
    %cfg.scalpsmooth = 2;
    cfg.brainthreshold = 0.25;
    cfg.scalpthreshold = 0.01;
    mri_seg_fem    = ft_volumesegment(cfg, mri_spm);
    
    save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_fem'],'mri_seg_fem');
    
end