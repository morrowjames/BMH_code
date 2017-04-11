% clear;close all; clc;

%ID = {'004'; '005'; '007'; '008'; '009'; '012'; '013'; '014'; '015'; '016'; '017'; '019'; '020'; '021'; '022';}; %'001';'002'
ID = {'004'}; %'001';'002'

%Need to re-do = {'001'}  

%Succesful re-do = 005 019
%pathIn = '/Volumes/OUTPUT DATA/Sian/EEG_coords_reconstruction/';
%pathOut = '/Volumes/OUTPUT DATA/Sian/EEG_coords_reconstruction/';
pathIn = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';
pathOut = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';

for a = 1:size(ID,1)

    load([pathIn,ID{a,1},filesep,ID{a,1},'_mri_acpc_reslice_segment_fem_mesh.mat']);

    
    %% Visualize surfaces
    figure
    ft_plot_mesh(mesh_fem(1),'facecolor','none','surfaceonly','true') % grey
    figure
    ft_plot_mesh(mesh_fem(2),'facecolor','none') % white
    figure
    ft_plot_mesh(mesh_fem(3),'facecolor','none') % csf
    figure
    ft_plot_mesh(mesh_fem(4),'facecolor','none') % skull
    figure
    ft_plot_mesh(mesh_fem(5),'facecolor','none') % scalp
    
    %% Visualize full mesh
    figure
    ft_plot_mesh(mesh_fem(3), 'facecolor',[0.2 0.2 0.2], 'facealpha', 0.3, 'edgecolor', [1 1 1], 'edgealpha', 0.05);
    hold on;
    ft_plot_mesh(mesh_fem(2),'edgecolor','none','facealpha',0.4);
    hold on;
    ft_plot_mesh(mesh_fem(1),'edgecolor','none','facecolor',[0.4 0.6 0.4]);

 
end