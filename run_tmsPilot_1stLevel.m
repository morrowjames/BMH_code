clear; clc; close all;

%addpath('/gpfs/M2Home/projects/Monash076/Nigel/TMS-fMRI/scripts/');

id = {'SC'};
cond = {'suprathreshold';'NOCOIL'; 'Coil'};
N = [186,186,186];
pathName = '/gpfs/M2Home/projects/Monash076/Morrowj/TMS-fMRI/';

for x = 1:size(id,1)
    tmsPilot_1stLevel(id{x,1},cond,N,pathName);
end