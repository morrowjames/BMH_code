function [] = tmsPilot_1stLevel(id,cond,N,pathName)

% This script performs the following steps
% 1. fMRI model specification
% 2. Model estimation
% 3. Contrast manager

% Inputs

% id = 'string'; Participant ID. Example: 'JC'.
% cond = {'cell'}; TMS-fMRI conditions. Example: {'LI1';'HI1'};
% N = [vector]; Number of acquisitions. Example: [60,66,66];
% pathName = 'string'; Path containing partipant ID folder;

% Dependencies

% The script assumes five low and five high intensity TMS conditions that have
% undergone preprocessing. Both folders are required 
% in a folder named as the participant [pilot,id]. A 'Timing' folder containing .m files with the
% onsets and durations is also required. The folders should be within the
% following path: /gpfs/M2Home/projects/Monash076/Nigel/TMS-fMRI/[pilot,id]/

% Created by Nigel Rogasch, 3-11-2016

spm_jobman('initcfg')

% Check if OUTPUT folder exists
% if ~isdir([pathName,id,'/OUTPUT']);
%     mkdir([pathName,id,'/OUTPUT']);
% end

%-----------------------------------------------------------------------
% Job saved on 02-Nov-2016 14:46:14 by cfg_util (rev $Rev: 6134 $)
% spm SPM - SPM12 (6225)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.stats.fmri_spec.dir = {[pathName,id,'/OUTPUT']};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.76;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 44;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
%%

%Load data
for i = 1:length(cond)
    for x = 1:N(1,i)
        matlabbatch{1}.spm.stats.fmri_spec.sess(i).scans{x,1} = [pathName,id,'/',id,'_',cond{i,1},'/swa',id,'_',cond{i,1},'.nii,',num2str(x)];
    end

    matlabbatch{1}.spm.stats.fmri_spec.sess(i).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi = {[pathName,id,'/Timing/stimuli_',id,'_',cond{i,1},'.mat']};
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).regress = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).multi_reg = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess(i).hpf = 128;
end
%%

% fMRI model specification
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

% Model Estimation 
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Contrast Manager
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'LI grt HI';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'HI grt LI';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 0;

spm_jobman('run',matlabbatch);
