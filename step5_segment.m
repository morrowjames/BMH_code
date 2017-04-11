
clear;close all; clc;

ID = {'002'}; %'001';'002'

pathIn = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
pathOut = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';

for a = 1:size(ID,1)

    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice']);
    
    %% Segment the anatomical data


    cfg           = [];
    cfg.output    = {'brain' 'skull' 'scalp'};
    cfg.brainsmooth = 5; %5
    %cfg.scalpsmooth = 2;
    cfg.brainthreshold = 0.4; %0.4
    cfg.scalpthreshold = 0.02; %0.01
    mri_seg_bem    = ft_volumesegment(cfg, mri_spm);

    %% Check segmentation brain
    cfg = [];
    cfg.method = 'ortho';
    cfg.funparameter = 'brain'; 
    ft_sourceplot(cfg,mri_seg_bem)

    %% Check segmentation skull
    cfg = [];
    cfg.method = 'ortho';
    cfg.funparameter = 'skull'; 
    ft_sourceplot(cfg,mri_seg_bem)

    %% Check segmentation scalp
    cfg = [];
    cfg.method = 'ortho';
    cfg.funparameter = 'scalp'; 
    ft_sourceplot(cfg,mri_seg_bem)
    
    save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment'],'mri_seg_bem');
end