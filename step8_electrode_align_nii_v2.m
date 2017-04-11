clear;close all; clc;

%ID = {'004'; '005'; '007'; '008'; '009'; '012'; '013'; '014'; '015'; '016'; '017'; '019'; '020'; '021'; '022';}; %'001';'002'
ID = {'005'};

pathIn = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
pathOut = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
% pathIn = '/Volumes/NO NAME/';
% pathOut = '/Volumes/NO NAME/';

for a = 1:size(ID,1)

    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_mesh_bem']);

    %% (1) Read electrode position from Localizer file

    % elec = ft_read_sens(experiment.electrodes);

    % elec = ft_read_sens('D:\nmda_tms_eeg\NEURONAV\Nigel_Study\DeDe\Debora_19900101_DeDo_4d87d5c3\Sessions\Session_20150922104459961\EEG\EEGMarkers20150922155315781.xml');

    elec.label = importElecNames([pathOut,filesep,ID{a,1},filesep,ID{a,1},'_electrode_positions.xlsx'],'Sheet1','A2:A66');
    elec.elecpos = importElecPos([pathOut,filesep,ID{a,1},filesep,ID{a,1},'_electrode_positions.xlsx']);
    elec.chanpos = importElecPos([pathOut,filesep,ID{a,1},filesep,ID{a,1},'_electrode_positions.xlsx']);
    
    %% (3) If the electrodes are not aligned to the head surface, interactively align them 
% (according to our visual judgement...)

    cfg = [];
    cfg.method = 'interactive';
    cfg.elec = elec;

    cfg.headshape = headmodel_bem.bnd(3); % scalp
    %cfg.headshape.pnt = cfg.headshape.pnt4;
    elec_aligned = ft_electroderealign(cfg);
    
    
    figure
    %plot scalp
    ft_plot_mesh(headmodel_bem.bnd(3),'edgecolor','none','facealpha',0.8,'facecolor',[0.6 0.6 0.8])
    hold on
    % plot electrodes
    ft_plot_sens(elec_aligned,'style','sk')
    
    save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_mesh_bem_elecAlign'],'elec_aligned');
    
end