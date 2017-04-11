
clear;close all; clc;

%ID = {'004'; '005'; '007'; '008'; '009'; '012'; '013'; '014'; '015'; '016'; '017'; '019'; '020'; '021'; '022';}; %'001';'002'
ID = {'003'}; %'001';'002'

%Need to re-do = {'001'}  

%Succesful re-do = 005 019
%pathIn = '/Volumes/OUTPUT DATA/Sian/EEG_coords_reconstruction/';
%pathOut = '/Volumes/OUTPUT DATA/Sian/EEG_coords_reconstruction/';
pathIn = '/Volumes/OUTPUT DATA/Sian/EEG_coords_reconstruction/';
pathOut = '/Volumes/OUTPUT DATA/Sian/EEG_coords_reconstruction/';

for a = 1:size(ID,1)

    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment']);

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
    
    save([pathOut,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_mesh'],'mesh_bem');
 
end
