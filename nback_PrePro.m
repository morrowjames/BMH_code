function [] = nback_PrePro(id,timePoint,sessType,N)

% This script performs the following steps
% 1. Slice timing correction
% 2. Realignment of images
% 3. Coregistration of EPI to T1
% 4. Normalisation of EPI to MNI
% 5. Smoothing

% Inputs

% id = 'string'; Participant ID. Example: 'P003'.
% timePoint = 'string'; 'PRE' | 'POST'.
% sessType = 'string'; Name of file type. Example: 'NBACK1'.
% N = integer; Number of acquisitions. Example: 119

% Dependencies

% The script assumes two nback sessions included in a folder titled 'nback'
% and a t1 image in a folder titled 't1'. Both folders are required 
% in a folder named as the participant id. The folder should be within the
% following path: /gpfs/M2Home/projects/Monash076/Nigel/gaba-tbs/

% Created by Nigel Rogasch, 26-10-2016

spm_jobman('initcfg')

%-----------------------------------------------------------------------
% Job saved on 31-Oct-2016 11:08:39 by cfg_util (rev $Rev: 6134 $)
% spm SPM - SPM12 (6225)
%-----------------------------------------------------------------------
%%

%Write out nii file names
for x = 1:N
    matlabbatch{1}.spm.temporal.st.scans{1,1}{x,1} = ['/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/',id,'/',timePoint,'_',sessType,'/',timePoint,'_',sessType,'.nii,',num2str(x)];
end

%%

% 1. Slice timing correction
matlabbatch{1}.spm.temporal.st.nslices = 44;
matlabbatch{1}.spm.temporal.st.tr = 2.50;
matlabbatch{1}.spm.temporal.st.ta = 2.443181; %TA = ((Slices-1/ Slices)*TR)
matlabbatch{1}.spm.temporal.st.so = [2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 42 44 1 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39 41 43];
matlabbatch{1}.spm.temporal.st.refslice = 1;
matlabbatch{1}.spm.temporal.st.prefix = 'a';

%2. Realingment of images
matlabbatch{2}.spm.spatial.realign.estwrite.data{1}(1) = cfg_dep('Slice Timing: Slice Timing Corr. Images (Sess 1)', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.sep = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.interp = 2;
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.eoptions.weight = '';
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.which = [2 1];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.interp = 4;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.mask = 1;
matlabbatch{2}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

% 3. Coregistration of EPI to T1
t1Ref = ['/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/',id,'/',timePoint,'_T1/',timePoint,'_T1.nii,1'];
matlabbatch{3}.spm.spatial.coreg.estwrite.ref = {t1Ref};
matlabbatch{3}.spm.spatial.coreg.estwrite.source(1) = cfg_dep('Realign: Estimate & Reslice: Mean Image', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','rmean'));
matlabbatch{3}.spm.spatial.coreg.estwrite.other(1) = cfg_dep('Realign: Estimate & Reslice: Realigned Images (Sess 1)', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','sess', '()',{1}, '.','cfiles'));
matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{3}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.interp = 4;
matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{3}.spm.spatial.coreg.estwrite.roptions.prefix = 'c';

% 4. Normalisation of EPI to MNI
matlabbatch{4}.spm.spatial.normalise.estwrite.subj(1).vol = {t1Ref};
matlabbatch{4}.spm.spatial.normalise.estwrite.subj(1).resample(1) = cfg_dep('Coregister: Estimate & Reslice: Coregistered Images', substruct('.','val', '{}',{3}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','cfiles'));
matlabbatch{4}.spm.spatial.normalise.estwrite.subj(2).vol = {t1Ref};
matlabbatch{4}.spm.spatial.normalise.estwrite.subj(2).resample = {t1Ref};
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.tpm = {'/usr/local/spm12/matlab2014a.r5556/tpm/TPM.nii'};
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{4}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                             78 76 85];
matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
matlabbatch{4}.spm.spatial.normalise.estwrite.woptions.interp = 4;

% 5. Smoothing
matlabbatch{5}.spm.spatial.smooth.data(1) = cfg_dep('Normalise: Estimate & Write: Normalised Images (Subj 1)', substruct('.','val', '{}',{4}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('()',{1}, '.','files'));
matlabbatch{5}.spm.spatial.smooth.fwhm = [6 6 6];
matlabbatch{5}.spm.spatial.smooth.dtype = 0;
matlabbatch{5}.spm.spatial.smooth.im = 0;
matlabbatch{5}.spm.spatial.smooth.prefix = 's';

spm_jobman('run',matlabbatch);

fprintf('%s %s %s is finished\n',id,timePoint,sessType);
