clear; close all; clc;

addpath ('/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/fieldtrip/fieldtrip-20150828');
ft_defaults;

% ID = {'H1';'H2'; 'H3'; 'H4'; 'H5'; 'H6'; 'H8'; 'H9'; 'H10'; 'H11'; 'H12';'H13'; 'H14'; 'H15'; 'H16'; 'H18'; 'H19'; 'H20'; 'H21';'H32'};
ID = {'004'};

% group = 'sz';

pathin = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/input/EEG_cleaned/';
pathout = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/output/';

site = {'PFC';'PAR'};
cond = {'H';'NH'};

for a=1:size(ID,1)
    for b = 1:size(site,1)
        for c = 1:size(cond,1)
            ERP = [];
            path = [pathin,ID{a,1},filesep];
            filename = [ID{a,1},'_',site{b,1},'_',cond{c,1},'_ep_ds_mark_clean_ica1_clean_ica2_clean_avref.set'];
            EEG = pop_loadset('filename', filename, 'filepath', path);
            ERP = eeglab2fieldtrip(EEG,'timelockanalysis');
            ERP.dimord = 'chan_time';
            savefile = [pathout,ID{a,1},'_ERP_',site{b,1},'_',cond{b,1},'.mat'];
            save (savefile, 'ERP');
        end
    end
end

