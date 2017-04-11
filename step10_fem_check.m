clear;close all; clc;

%ID = {'004'; '005'; '007'; '008'; '009'; '012'; '013'; '014'; '015'; '016'; '017'; '019'; '020'; '021'; '022';}; %'001';'002'
ID = {'012'};

pathIn = '/Volumes/OUTPUT DATA/Sian/EEG_coords_reconstruction/';
pathOut = '/Volumes/OUTPUT DATA/Sian/EEG_coords_reconstruction/';
% pathIn = '/Volumes/NO NAME/';
% pathOut = '/Volumes/NO NAME/';

for a = 1:size(ID,1)
    
    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_fem.mat']);

    segmentation = ft_datatype_segmentation(mri_seg_fem,'segmentationstyle','indexed');

    cfg              = [];
    cfg.funparameter = 'seg';
    cfg.funcolormap  = lines(6); % distinct color per tissue
    cfg.location     = 'center';
    cfg.atlas        = segmentation;    % the segmentation can also be used as atlas 
    ft_sourceplot(cfg, segmentation);

end