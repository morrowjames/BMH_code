clear;close all; clc;

ID = {'002';}; %DONT RUN 010 011 018!!!


addpath ('/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/fieldtrip/fieldtrip-20150828');
ft_defaults;

for a = 1:size(ID,1)
    pathName = ['/gpfs/M2Home/projects/Monash076/Morrowj/DICOM/',ID{a,1},'_Dicom/'];
    pathOut = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output';

    %% Choose subject's MRI data

    pathIn = ['/gpfs/M2Home/projects/Monash076/Morrowj/DICOM/',ID{a,1},'_Dicom/'];
    nifti = '002.nii';
    experiment = [pathIn,nifti]; % DICOM or Nifti dataset

    %% Read  MRI data
    mri = ft_read_mri(experiment); % Dicom: reads all slides automatically (give just 1st one, is enough)
    
    mkdir(pathOut,ID{a,1});
   
    save([pathOut,ID{a,1},filesep,ID{a,1},'_raw_mri'],'mri');
end