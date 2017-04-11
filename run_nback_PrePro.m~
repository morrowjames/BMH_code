clear; clc; close all;

addpath('/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/scripts');

id = {'coil_test'};
%timePoint = {'PRE';'POST'};
timePoint = {'32'; '8'};
sessType = {'REST'};
% sessType = {'NBACK2'};
N = 186;

for x = 1:size(id,1)
    for y = 1:size(timePoint,1)
        for z = 1:size(sessType,1)
            nback_PrePro(id{x,1},timePoint{y,1},sessType{z,1},N);
        end
    end
end