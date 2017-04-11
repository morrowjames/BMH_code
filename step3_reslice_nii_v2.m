
clear;close all; clc;

ID = {'005'}; %'001';'002'

pathIn = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
pathOut = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';


for a = 1:size(ID,1)
    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc']);
    
    %% Determine the handedness of the transformation matrix

    %  If the transformation matrix is describing the transformation from voxel space to a world coordinate system 
    %  where the XYZ-axes are described with right-handed coordinate axes, AND the transformation matrix is also right-handed, 
    %  then 'right is right, and left is left'.
    % 
    %  If the transformation matrix is describing the transformation from voxel space to a world coordinate system
    %  where the XYZ-axes are described with right-handed coordinate axes, AND the transformation matrix is left-handed, 
    %  then 'right is left, and left is right', and equivalently 'the image is flipped'.

    % If handedness has a positive value, the transformation-matrix is right-handed, otherwise it's left-handed. 


    handedness = det(mri_spm.transform(1:3,1:3))

    %% Re-slice the MRI volume so that all the 3 dimensions have the same number of voxels
    % (it would be also possible to specify  cfg.dim = [nx ny nz]; )

    cfg = [];  
    cfg.resolution = 1; % 1 mm thick slices
    cfg.dim = [256 256 256];
    % cfg.dim = [300 300 300];

    mri_spm = ft_volumereslice(cfg,mri_spm);
    
    save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice'],'mri_spm');
end