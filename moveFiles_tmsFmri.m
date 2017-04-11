clear; clc; close all;

pathName = '/gpfs/M2Home/projects/Monash076/Morrowj/TMS-fMRI/';

ID = {'SC'};

darisID = '1008.2.72';

newFolders = {'SC_T1';'SC_suprathreshold';'SC_subthreshold';'SC_subthresholdmove';'SC_NOCOIL'; 'SC_Coil'};

for i = 1:size(ID,1)
    % Make new directories
    for x = 1:size(newFolders)
        mkdir([pathName,ID{i,1}],newFolders{x,1});
    end

    % Find folders with processed files
    temp = dir([pathName,ID{i,1},filesep,darisID]);
    dirNames = {temp.name};
    [~,ind] = find(strncmp(darisID,dirNames,size(darisID,2)));
    num = dirNames{1,ind}(1,size(darisID,2)+2:end);
    filePathNum = [darisID,'.',num];

    pathNameFull = [pathName,ID{i,1},filesep,darisID,filesep,[darisID,'.',num],filesep,[darisID,'.',num,'.1'],filesep,[darisID,'.',num,'.1.4']];
    allFiles = dir(pathNameFull);
    allFilesNames = {allFiles.name}';

    for x = 1:size(newFolders,1)
        for y = 3:size(allFiles,1)
            tempFiles = dir([pathNameFull,filesep,allFilesNames{y,1}]);
            tempFilesNames = {tempFiles.name}';
            if sum(strcmp([newFolders{x,1},'.nii'],tempFilesNames)) > 0;
                copyfile([pathNameFull,filesep,allFilesNames{y,1},filesep,newFolders{x,1},'.nii'],[pathName,ID{i,1},filesep,newFolders{x,1}]);
            end
        end
    end
end