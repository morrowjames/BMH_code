%%
%clear;close all; clc;

%ID = {'004'; '005'; '007'; '008'; '009'; '012'; '013'; '014'; '015'; '016'; '017'; '019'; '020'; '021'; '022';}; %'001';'002'
ID = {'004'};

pathIn = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
pathOut = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
% pathIn = '/Volumes/NO NAME/';
% pathOut = '/Volumes/NO NAME/';
%%
for a = 1:size(ID,1)
   %% Prepare leadfield FEM headmodel
load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_fem_mesh_headmodel.mat']);
load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_mesh_bem_elecAlign.mat']);
        cfg = [];
        cfg.elec = elec_aligned;
        cfg.grid = sourcemodel;
        cfg.headmodel = headmodel_fem;
        cfg.normalize = 'yes'; % very important!!! (use default vaue 0.5)
        leadfield_fem = ft_prepare_leadfield(cfg);
    
    save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_fem_mesh_headmodel_leadfield'],'headmodel_fem');

end