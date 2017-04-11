clear; close all; clc

pathIn = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
pathOut = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';

% NOTE: the path to the template file is user-specific
load('/gpfs/M2Home/projects/Monash076/Morrowj/DICOM/002_Dicom/002.nii')
ftpath   = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/fieldtrip/fieldtrip-20150828'; % this is the path to fieldtrip at Donders


template = ft_read_mri('002.nii');
template.coordsys = 'spm'; % so that FieldTrip knows how to interpret the coordinate system

cfg = [];
cfg.method = 'ortho';
ft_sourceplot(cfg,template)


cfg = [];  
cfg.resolution = 1; % 1 mm thick slices
cfg.dim = [256 256 256];

template = ft_volumereslice(cfg,template);

 
% segment the template brain and construct a volume conduction model (i.e. head model): 
% this is needed to describe the boundary that define which dipole locations are 'inside' the brain.
cfg          = [];
cfg.output    = {'brain' 'skull' 'scalp'};
template_seg = ft_volumesegment(cfg, template);


%% Check segmentation brain
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'brain'; 
ft_sourceplot(cfg,template_seg)

%% Check segmentation skull
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'skull'; 
ft_sourceplot(cfg,template_seg)

%% Check segmentation scalp
cfg = [];
cfg.method = 'ortho';
cfg.funparameter = 'scalp'; 
ft_sourceplot(cfg,template_seg)


%%
 
cfg          = [];
cfg.method   = 'singleshell';
template_headmodel = ft_prepare_headmodel(cfg, template_seg);
%template_headmodel = ft_convert_units(template_headmodel, 'cm');

 
% construct the dipole grid in the template brain coordinates
% the negative inwardshift means an outward shift of the brain surface for inside/outside detection
cfg = [];
cfg.grid.xgrid  = -250:6:250;
cfg.grid.ygrid  = -250:6:250;
cfg.grid.zgrid = -250:6:250;
cfg.headmodel   = template_headmodel;
cfg.elec = hdr.elec;
cfg.inwardshift = -25;
cfg.grid.unit = 'mm';
cfg.grid.tight = 'yes';
template_grid = ft_prepare_sourcemodel(cfg);
 
% make a figure with the template head model and dipole grid
figure
hold on
ft_plot_vol(template_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
ft_plot_mesh(template_grid.pos(template_grid.inside,:));
ft_plot_sens(hdr.elec,'label','label','style','m*','Markersize',10);

save template_grid template_grid