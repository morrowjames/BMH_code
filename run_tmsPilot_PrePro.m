clear; clc; close all;

addpath('/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/scripts/');

id = {'SC'};
condition = {'suprathreshold';'NOCOIL'; 'Coil'};
N = [186, 186, 186]; % Number of volumes/acquistions
TR = [2.76,2.76,2.76]; % TR
TA = [2.5995,2.5995,2.5995]; % TA = TR-(TR/Nslices); Note that with delay TR should be acquistion time
pathName = '/gpfs/M2Home/projects/Monash076/Morrowj/TMS-fMRI/';

% id = {'JC3'};
% condition = {'HI4';'HI5';'LIMV1';'LIMV2';'LIMV3'};
% N = [66,66,66,66,66]; % Number of volumes/acquistions
% TR = [2.76,2.76,2.76,2.76,2.76]; % TR
% TA = [2.5995,2.5995,2.5995,2.5995,2.5995]; % TA = TR-(TR/Nslices); Note that with delay TR should be acquistion time

for x = 1:size(id,1)
    for y = 1:size(condition,1)
        tmsPilot_PrePro(id{x,1},condition{y,1},N(1,y),TR(1,y),TA(1,y),pathName);
    end
end