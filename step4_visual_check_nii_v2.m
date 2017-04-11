clear;close all; clc;

ID = {'005'}; %'001';'002'

pathIn = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
pathOut = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';

for a = 1:size(ID,1)
    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice']);

    %% Visualize MRI data to check that everything is ok

    cfg = [];
    cfg.method = 'ortho';
    ft_sourceplot(cfg,mri_spm)
    
end