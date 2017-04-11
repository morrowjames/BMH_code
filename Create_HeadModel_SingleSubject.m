


%%% ----------------  CREATE HEAD MODEL FOR SOURCE LOCALIZATION TEPs ------------------- %%%%

% This process needs to be done by hand and  step by step for each subject

% The MRI data needs to be created just once 
% The elctrode position needs to be loaded and aligned to the MRI for each
% session separately


%% Choose subject's MRI data

experiment = EXPERIMENTS(27) % DICOM or Nifti dataset

%% Read  MRI data
mri = ft_read_mri(experiment.mri); % Dicom: reads all slides automatically (give just 1st one, is enough)


%% Visualize MRI data to check that everything is ok

cfg = [];
cfg.method = 'ortho';
ft_sourceplot(cfg,mri)


%% Convert anatomical MRI to MNI coordinate system (= RAS)

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'spm'; % MNI coordinate system (=RAS)
cfg.viewmode = 'ortho';
mri_spm = ft_volumerealign(cfg,mri);

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

%% Visualize MRI data to check that everything is ok

cfg = [];
cfg.method = 'ortho';
ft_sourceplot(cfg,mri_spm)



%% Re-slice the MRI volume so that all the 3 dimensions have the same number of voxels
% (it would be also possible to specify  cfg.dim = [nx ny nz]; )

cfg = [];  
cfg.resolution = 1; % 1 mm thick slices
cfg.dim = [256 256 256];

mri_spm = ft_volumereslice(cfg,mri_spm);

%% Visualize MRI data to check that everything is ok

cfg = [];
cfg.method = 'ortho';
ft_sourceplot(cfg,mri_spm)



%% Save mri data
save([experiment.path '\ForwardModel'],'mri_spm','-v7.3')

%save([experiment.path '\ForwardModel'],'mri_spm','-append')



%% Segment the anatomical data


cfg           = [];
cfg.output    = {'brain' 'skull' 'scalp'};
cfg.brainsmooth = 5;
%cfg.scalpsmooth = 2;
cfg.brainthreshold = 0.4;
cfg.scalpthreshold = 0.01;
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

%% If everything looks fine: Save segmented MRI

save([experiment.path '\ForwardModel'],'mri_seg_bem','-append')


%% Prepare mesh: creates surfaces at the boarders of the different tissue-types 
cfg = [];
cfg.method = 'projectmesh';
cfg.tissue = {'brain' 'skull' 'scalp'};
cfg.numvertices = [6000 4000 2000];
mesh_bem = ft_prepare_mesh(cfg, mri_seg_bem);

%% Visualize surfaces
figure
ft_plot_mesh(mesh_bem(1),'facecolor','none') % brain
figure
ft_plot_mesh(mesh_bem(2),'facecolor','none') % skull
figure
ft_plot_mesh(mesh_bem(3),'facecolor','none') % scalp

%% Visualize full mesh
figure
ft_plot_mesh(mesh_bem(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
hold on;
ft_plot_mesh(mesh_bem(2),'edgecolor','none','facealpha',0.4);
hold on;
ft_plot_mesh(mesh_bem(1),'edgecolor','none','facecolor',[0.4 0.6 0.4]);

%% Save mech for BEM model

save([experiment.path '\ForwardModel'],'mesh_bem','-append')

%% Prepare head model (volume conduction model)
cfg = [];
cfg.method = 'bemcp';
headmodel_bem = ft_prepare_headmodel(cfg,mesh_bem);                 

%% Save the head model (it contains already also the mesh in field bnd)

save([experiment.path '\ForwardModel'],'headmodel_bem','-append')

%% %%% -------------------------------- Align the EEG electrodes!!! -------------------------- %%%

%% (1) Read electrode position from Localizer file

elec = ft_read_sens(experiment.electrodes);



%% (2) Check if alignment is necessary
load([experiment.path '\ForwardModel'],'headmodel_bem');
figure
%plot scalp
ft_plot_mesh(headmodel_bem.bnd(3),'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8])
hold on
% plot electrodes
ft_plot_sens(elec,'style','sk')


%% (3) If the electrodes are not aligned to the head surface, interactively align them 
% (according to our visual judgement...)

cfg = [];
cfg.method = 'interactive';
cfg.elec = elec;

cfg.headshape = headmodel_bem.bnd(3); % scalp
elec_aligned = ft_electroderealign(cfg);



%%  Load hdr and order elec_aligned accordign to original order of signal acquisition

load([experiment.path '\Data'],'hdr')

 [~, original_order] = ismember(hdr.elec.label,elec_aligned.label);
 
 
elec_aligned.chanpos = elec_aligned.chanpos(original_order,:);
elec_aligned.elecpos = elec_aligned.elecpos(original_order,:);
elec_aligned.label = elec_aligned.label(original_order);

%% Check if the alignment is ok

figure
%plot scalp
ft_plot_mesh(headmodel_bem.bnd(3),'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8])
hold on
% plot electrodes
ft_plot_sens(elec_aligned,'style','sk')

%% Save electrodes aligned

save([experiment.path '\ForwardModel'],'elec_aligned','-append')


%% -------------------------------- Source model ---------------------------- %%

%% Make the individual subject's grid

load([experiment.path '\ForwardModel'],'elec_aligned','mri_spm','headmodel_bem')
load('template_grid')


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

save([experiment.path '\ForwardModel'],'sourcemodel','-append')




%% --------------------------  FORWARD SOLUTION  BEM------------------------------ %%

%% Prepare leadfield BEM
cfg = [];
cfg.elec = elec_aligned;
cfg.grid = sourcemodel;
cfg.headmodel = headmodel_bem;
leadfield_bem = ft_prepare_leadfield(cfg);

%% Save subject's forward solution (leadfield)

save([experiment.path '\ForwardModel'],'leadfield_bem','-append')

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

save([experiment.path '\ForwardModel'],'mri_seg_fem','-append')

%% Prepare mesh for FEM head model

cfg        = [];
cfg.shift  = 0.3;
cfg.method = 'hexahedral';
mesh_fem = ft_prepare_mesh(cfg,mri_seg_fem);

%% Save mesh for FEM model

save([experiment.path '\ForwardModel'],'mesh_fem','-append')

%% Prepare head model (volume conduction model)--> NB double check order of mesh.tissuelabel
load([experiment.path '\ForwardModel'],'mesh_fem')

cfg = [];
cfg.method ='simbio';
cfg.conductivity = [1.79 0.33 0.43 0.01 0.14];   % order follows mesh.tissuelabel: csf gray scalp skull white
headmodel_fem = ft_prepare_headmodel(cfg,mesh_fem);

%% Save the head model (it contains already also the mesh in field bnd)

save([experiment.path '\ForwardModel'],'headmodel_fem','-append')




%% Prepare leadfield FEM headmodel
load([experiment.path '\ForwardModel'],'sourcemodel','headmodel_fem','elec_aligned')

cfg = [];
cfg.elec = elec_aligned;
cfg.grid = sourcemodel;
cfg.headmodel = headmodel_fem;
cfg.normalize = 'yes'; % very important!!! (use default vaue 0.5)
leadfield_fem = ft_prepare_leadfield(cfg);

%% Save subject's forward solution (leadfield)

save([experiment.path '\ForwardModel'],'leadfield_fem','-append')


    
   















































