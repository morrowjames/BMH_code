clear;close all; clc;

ID = {'001'}; %'001';'002'

pathIn = '/Users/jamesmorrow/Desktop/Output/';
pathOut = '/Users/jamesmorrow/Desktop/Output/';

for a = 1:size(ID,1)
    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice']);

    %% Visualize MRI data to check that everything is ok

    cfg = [];
    cfg.method = 'ortho';
    ft_sourceplot(cfg,mri_spm)
    
end