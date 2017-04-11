clc clear all

% Add the path where the trigger scripts are
addpath = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/scripts/TMS_trigger';

% Open serial/com port connection. Default is COM1. Open mp_open_serial to
% change com port if needed. Once a port is opened, if this function is run
% again it will error
s = mp_open_serial;

% Read data from stimulator to make sure connection is opened. Matrix
% should start with 254 and end with 255 (or similar). Value in the middle
% will be stimulation intensity

msg_from_stimulator=fread(s);

%%

% commands to adjust and trigger stimualtor

mp_enable_stimulator(s,0); % 0 = disable, 1 = enable

mp_set_amplitude(s,30); % value = intensity

mp_single_trigger(s); % triggers one pulse. No input arguments except serial port

mp_start_train; % if pulse train is defined this will start

%%

%close serial port connection and delete port object. Solves the problem of
%trying to open a second port if you happen to do that

fclose(s);delete(s);clear s;






