clear; clc; close all;

% addpath('/gpfs/M2Home/projects/Monash076/Nigel/gaba-tbs/scripts/');
% 
% id = {'P002';'P003'};
% timePoint = {'PRE';'POST'};
% N = 119;

addpath('/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/scripts');

id = {'JM'};
%timePoint = {'PRE';'POST'};
timePoint = {'PRE'};
N = 119;

for x = 1:size(id,1)
    for y = 1:size(timePoint,1)
        nback_1stLevel(id{x,1},timePoint{y,1},N);
    end
end