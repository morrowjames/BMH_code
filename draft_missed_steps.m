%% -------------------------------- Source model ---------------------------- %%
% clear;close all; clc;

%ID = {'004'; '005'; '007'; '008'; '009'; '012'; '013'; '014'; '015'; '016'; '017'; '019'; '020'; '021'; '022';}; %'001';'002'
ID = {'004'};

pathIn = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
pathOut = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
% pathIn = '/Volumes/NO NAME/';
% pathOut = '/Volumes/NO NAME/';


a = 1:size(ID,1);
%% Make the individual subject's grid

    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_mesh_bem_elecAlign.mat']);

    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice.mat']);

    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_mesh_bem.mat']);

    load('template_grid');


cfg = [];
cfg.grid.warpmni = 'yes';
cfg.grid.template = template_grid; 
cfg.grid.nonlinear = 'yes'; % use non linear normalization
cfg.mri = mri_spm; % use subject's mri in mni coordinates
cfg.grid.unit = 'mm';
sourcemodel = ft_prepare_sourcemodel(cfg);



%% Plot the single subject'e head model and grid positions

figure
hold on
ft_plot_vol(headmodel_bem,'edgecolor','none','facealpha',0.3);
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));



%% Double check that source model is aligned to the segmented volume and to the electrodes
figure
hold on
%ft_plot_mesh(sourcemodel)
plot3(sourcemodel.pos(sourcemodel.inside,1),sourcemodel.pos(sourcemodel.inside,2),sourcemodel.pos(sourcemodel.inside,3),'.');
ft_plot_vol(headmodel_bem,'facecolor','brain','edgecolor','none','facealpha',0.3);
ft_plot_sens(elec_aligned,'label','label','style','y*','Markersize',10);

%%

headmodel_brain = headmodel_bem;
headmodel_brain = rmfield(headmodel_brain,'bnd');
headmodel_brain.bnd = headmodel_bem.bnd(1);

figure
hold on
ft_plot_vol(headmodel_brain, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:));
ft_plot_sens(elec_aligned,'label','label','style','m*','Markersize',10);

%% Double check headmodel for inside/outside brain
figure
hold on
ft_plot_vol(headmodel_grid,'facecolor','brain','edgecolor','none','facealpha',0.3);
ft_plot_vol(headmodel_brain, 'facecolor', 'skin', 'edgecolor', 'none');alpha 0.1; camlight;



%% Save subject's grid (sourcemodel)

save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_missed_steps'],'headmodel_bem','sourecemodel');


%% --------------------------  FORWARD SOLUTION  BEM------------------------------ %%

%% Prepare leadfield BEM
cfg = [];
cfg.elec = elec_aligned;
cfg.grid = sourcemodel;
cfg.headmodel = headmodel_bem;
leadfield_bem = ft_prepare_leadfield(cfg);

%% Save subject's forward solution (leadfield)

save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_missed_steps_2'],'leadfield_bem', 'headmodel_bem');


%%  Before doing the FEM headmodel (very time consuming!!) check whether the cortex is fully covered!
%%% In case ther are some not covered areas (expecially in the center or if
%%% these uncovered areas are big), preprocess the MRI with SPM12 and
%%% repeat the BEM model creation procedure... double check again if the
%%% cortex has "holes" and so on...if the preprocessing in SPM does not make any difference then check that the grid is well positioned!


% Call the script for checking that the cortex is fully covered (uses fist clenching data)
SourceLocalization_CheckIfCortexIsFullyCovered

%% --------------------------------- FEM HEAD MODEL and FORWARD SOLUTION ---------------- %%%%%%%%%%%%%%

%% Segment the anatomical data for FEM headmodel
% load([experiment.path '\Data'],'mri_spm')

cfg           = [];
cfg.output    = {'gray','white','csf','skull','scalp'};
cfg.brainsmooth = 2;
%cfg.scalpsmooth = 2;
cfg.brainthreshold = 0.25;
cfg.scalpthreshold = 0.01;
mri_seg_fem    = ft_volumesegment(cfg, mri_spm);

%% Check segmentation FEM
segmentation = ft_datatype_segmentation(mri_seg_fem,'segmentationstyle','indexed');

cfg              = [];
cfg.funparameter = 'seg';
cfg.funcolormap  = lines(6); % distinct color per tissue
cfg.location     = 'center';
cfg.atlas        = segmentation;    % the segmentation can also be used as atlas 
ft_sourceplot(cfg, segmentation);

%% If everything looks fine: Save segmented MRI

save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_missed_steps_3'],'mri_seg_fem')

%% Prepare mesh for FEM head model

cfg        = [];
cfg.shift  = 0.3;
cfg.method = 'hexahedral';
mesh_fem = ft_prepare_mesh(cfg,mri_seg_fem);

%% Save mesh for FEM model

save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_missed_steps_4'],'mesh_fem');

%% Prepare head model (volume conduction model)--> NB double check order of mesh.tissuelabel
% load([experiment.path '\ForwardModel'],'mesh_fem')

cfg = [];
cfg.method ='simbio';
cfg.conductivity = [1.79 0.33 0.43 0.01 0.14];   % order follows mesh.tissuelabel: csf gray scalp skull white
headmodel_fem = ft_prepare_headmodel(cfg,mesh_fem);

%% Save the head model (it contains already also the mesh in field bnd)

save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_missed_steps_5'],'headmodel_fem')




%% Prepare leadfield FEM headmodel
% load([experiment.path '\ForwardModel'],'sourcemodel','headmodel_fem','elec_aligned')

cfg = [];
cfg.elec = elec_aligned; 
cfg.grid = sourcemodel;
cfg.headmodel = headmodel_fem;
cfg.normalize = 'yes'; % very important!!! (use default vaue 0.5)
leadfield_fem = ft_prepare_leadfield(cfg);

% Save subject's forward solution (leadfield)

% save([experiment.path '\ForwardModel'],'leadfield_fem','-append')
save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_missed_steps_leadfield'],'leadfield_fem')

