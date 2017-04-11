function [] = nback_1stLevel(id,timePoint,N)

% This script performs the following steps
% 1. fMRI model specification
% 2. Model estimation
% 3. Contrast manager

% Inputs

% id = 'string'; Participant ID. Example: 'P003'.
% timePoint = 'string'; 'PRE' | 'POST'.
% N = integer; Number of acquisitions. Example: 44

% Dependencies

% The script assumes two nback sessions included in a folder titled 'NBACK' that have
% undergone preprocessing. Both folders are required 
% in a folder named as the participant id. A 'Timing' folder containing .m files with the
% onsets and durations is also required. The folders should be within the
% following path: /gpfs/M2Home/projects/Monash076/Nigel/gaba-tbs/

% Created by Nigel Rogasch, 31-10-2016

spm_jobman('initcfg')

%-----------------------------------------------------------------------
% Job saved on 31-Oct-2016 15:42:53 by cfg_util (rev $Rev: 6134 $)
% spm SPM - SPM12 (6225)
%-----------------------------------------------------------------------

% Check if OUTPUT folder exists
if ~isdir(['/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/',id,'/OUTPUT_',timePoint]);
    mkdir(['/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/',id],['OUTPUT_',timePoint]);
end

% fMRI model specification - session 1
matlabbatch{1}.spm.stats.fmri_spec.dir = {['/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/',id,'/OUTPUT_',timePoint]};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2.5;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 44;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
%%
for x = 1:N
    matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans{x,1} = ['/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/',id,'/',timePoint,'_NBACK1/swa',timePoint,'_NBACK1.nii,',num2str(x)];
end

%%
matlabbatch{1}.spm.stats.fmri_spec.sess(1).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi = {['/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/',id,'/Timing/',id,'_',timePoint,'_NBACK1.mat']};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi_reg = {['/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/',id,'/',timePoint,'_NBACK1/rp_a',timePoint,'_NBACK1.txt']};
matlabbatch{1}.spm.stats.fmri_spec.sess(1).hpf = 128;
%%

% fMRI model specification - session 2
for x = 1:N
    matlabbatch{1}.spm.stats.fmri_spec.sess(2).scans{x,1} = ['/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/',id,'/',timePoint,'_NBACK2/swa',timePoint,'_NBACK2.nii,',num2str(x)];
end

%%
matlabbatch{1}.spm.stats.fmri_spec.sess(2).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi = {['/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/',id,'/Timing/',id,'_',timePoint,'_NBACK2.mat']};
matlabbatch{1}.spm.stats.fmri_spec.sess(2).regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess(2).multi_reg = {['/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/',id,'/',timePoint,'_NBACK2/rp_a',timePoint,'_NBACK2.txt']};
matlabbatch{1}.spm.stats.fmri_spec.sess(2).hpf = 128;

% fMRI model specification 
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{1}.spm.stats.fmri_spec.mask = {''};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'AR(1)';

% Model estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

% Contrast manager
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('fMRI model specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = '0back grt 2back';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'repl';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = '2back grt 0back';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [-1 1];
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'repl';
matlabbatch{3}.spm.stats.con.delete = 0;

spm_jobman('run',matlabbatch);
