clear;close all; clc;

ID = {'005'}; %'001';'002'

pathIn = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
pathOut = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';

for a = 1:size(ID,1)
    load([pathIn,ID{a,1},filesep,ID{a,1},'_raw_mri']);
    
    %% Convert anatomical MRI to MNI coordinate system (= RAS)

    cfg = [];
    cfg.method = 'interactive';
    cfg.coordsys = 'spm'; % MNI coordinate system (=RAS)
    cfg.viewmode = 'ortho';
    mri_spm = ft_volumerealign(cfg,mri);
    
    save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc'],'mri_spm');
end
