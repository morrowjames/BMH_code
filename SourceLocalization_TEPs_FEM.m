

%%  %%%----------------------------------- SOURCE TEPs EXCITABILITY SINGLE SUBJECT ------------------------- %%%%%


%%
clear; close all; clc;

addpath ('/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/fieldtrip/fieldtrip-20150828');
ft_defaults;

% ID = {'H1';'H2'; 'H3'; 'H4'; 'H5'; 'H6'; 'H8'; 'H9'; 'H10'; 'H11'; 'H12';'H13'; 'H14'; 'H15'; 'H16'; 'H18'; 'H19'; 'H20'; 'H21';'H32'};
ID = {'004'};

% group = 'sz';

pathin = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
pathout = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';



tic
for a=1:size(ID,1) %experiment = EXPERIMENTS([1 3 4 5 8 14 16 19 20 21 22 27])

     initialVars = who; % this is to keep the common variables in the workspace
%%% Load data and forward model variables
% warning(['Starting Subject ', experiment.path(13:16),'...'])
load([pathin,ID{a,1},filesep,ID{a,1},'_ERP_PFC_H.mat']);
load([pathout,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_mesh_bem_elecAlign.mat'])
load([pathout,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_missed_steps_leadfield.mat'])
load([pathout,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_fem_mesh_headmodel.mat'])

% load([experiment.path '\ForwardModel'],'elec_aligned','leadfield_fem','headmodel_fem','mri_spm')

% Remove last three electrodes
elec_aligned.label(63:65,:) = [];
elec_aligned.chanpos(63:65,:) = [];
elec_aligned.elecpos(63:65,:) = [];

%%% Source Analysis: Create spatial filter using the LCMV beamformer

cfg                  = [];
cfg.method           = 'lcmv';
cfg.elec             = elec_aligned;
cfg.grid             = leadfield_fem; % leadfield, which has the grid information
cfg.headmodel        = headmodel_fem; % volume conduction model (headmodel)
cfg.keepfilter       = 'yes';
SRC_ERP              = ft_sourceanalysis(cfg, ERP);
% SRC_peak_stim        = ft_sourceanalysis(cfg, ERP_peak_stim);
% SRC_peak_nostim      = ft_sourceanalysis(cfg, ERP_peak_nostim);
% SRC_trough_stim      = ft_sourceanalysis(cfg, ERP_trough_stim);
% SRC_trough_nostim    = ft_sourceanalysis(cfg, ERP_trough_nostim);

%%% Save source data

save([pathout,ID{a,1},filesep,ID{a,1},'_SRC_ERP']);

end
% %%
% 
% %%% Extract source timecourse with ft_sourcedescriptives
% 
% cfg = [];
% cfg.projectmom         = 'yes';
% cfg.fixedori           = 'over_trials';
% SRC_baselineTC         = ft_sourcedescriptives(cfg,SRC_baseline);
% SRC_peak_stimTC        = ft_sourcedescriptives(cfg, SRC_peak_stim);
% SRC_peak_nostimTC      = ft_sourcedescriptives(cfg, SRC_peak_nostim);
% SRC_trough_stimTC      =ft_sourcedescriptives(cfg, SRC_trough_stim);
% SRC_trough_nostimTC    = ft_sourcedescriptives(cfg, SRC_trough_nostim);
% 
% 
% Npos = size(SRC_baselineTC.pos,1);
% Ntime = length(SRC_baselineTC.time);
% 
% SRC_baselineTC.avg.momint = nan(Npos,Ntime);
% SRC_peak_stimTC.avg.momint = nan(Npos,Ntime);
% SRC_peak_nostimTC.avg.momint = nan(Npos,Ntime);
% SRC_trough_stimTC.avg.momint = nan(Npos,Ntime);
% SRC_trough_nostimTC.avg.momint = nan(Npos,Ntime);
% 
% idx = find(SRC_baselineTC.inside);
% 
% for i = 1:length(idx)
%     index = idx(i);
%     SRC_baselineTC.avg.momint(index,:) = SRC_baselineTC.avg.mom{index};
%     SRC_peak_stimTC.avg.momint(index,:) = SRC_peak_stimTC.avg.mom{index};
%     SRC_peak_nostimTC.avg.momint(index,:) = SRC_peak_nostimTC.avg.mom{index};
%     SRC_trough_stimTC.avg.momint(index,:) = SRC_trough_stimTC.avg.mom{index};
%     SRC_trough_nostimTC.avg.momint(index,:) = SRC_trough_nostimTC.avg.mom{index};
% end
% 
% %%% Save time course source data
% 
% save([experiment.path '\Data'],'SRC_baselineTC','SRC_peak_stimTC','SRC_peak_nostimTC','SRC_trough_stimTC','SRC_trough_nostimTC','-append')
% 
% 
%  clearvars('-except',initialVars{:})
% end
toc
% 
% 
% 
