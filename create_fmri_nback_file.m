clear; clc; close all;

ID = 'JM';

pathName = '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/';

cond = {'PRE_NBACK1';...
    'PRE_NBACK2'};...
    %'POST_NBACK1';...
    %'POST_NBACK2'};

fileName = {'NBACK_1-999-1.txt';...
    'NBACK_2-999-2.txt'};
   %'CogBlock1-003-2.txt';...
   %'CogBlock2-003-2.txt'};
    
for i = 1:size(cond,1)
    %Read in data
    [Titles,Header,Start,VarName4] = importfileNback([pathName,ID,filesep,'Timing',filesep,fileName{i,1}]);

    titleNames = {'ZeroBackSlide';'TwoBackSlide'};
    titleNamesShort = {'0back';'2back'};
    
    for z = 1:size(titleNames,1)
        %Find onset times and trial length times
        output = [];
        trialLog = zeros(size(Titles,1),1);
        for x = 1:size(Titles,1)
            output(x,1) = any(regexp(Titles{x,1},titleNames{z,1})) && any(regexp(Titles{x,1},'OnsetTime:$')) && ~any(regexp(Titles{x,1},'ToOnsetTime:$'));
        end

        clear zeroTime
        onsTime = Header(logical(output));
        temp = (1:1:size(Start,1))';
        onsInd = temp(logical(output));
        
%         [trialTime,~] = find(strncmp('TrialTime:',Titles));
        for x = 1:size(Titles,1)
            if isempty(regexp(Titles{x,1},regexptranslate('wildcard','*TrialTime:*')));
                TitCount(x,1) = 0;
            else
                TitCount(x,1) = 1;
            end
        end
        [trialTime,~] = find(TitCount == 1);

        for x = 1:size(onsTime,1)
            onsTimeSec(x,1) = str2num(onsTime{x,1})./1000;
            [~,trialInd] = min(abs(onsInd(x,1)-trialTime));
            trialTimeSec(x,1) = str2num(Header{trialTime(trialInd,1),1})./1000;
        end

        %Find start and end of each block
        diffOns = diff(onsTimeSec);
        ons = [];
        off = [];
        ons = round(onsTimeSec(1,1));
        for x = 1:size(diffOns,1)
            if diffOns(x,1) > 4
                ons(1,size(ons,2)+1) = onsTimeSec(x+1,1);
                off(1,size(off,2)+1) = onsTimeSec(x,1)+trialTimeSec(x,1);
            end
        end
        off(1,size(off,2)+1) = onsTimeSec(end,1)+trialTimeSec(end,1);
        
        %Calculation durations
        dur = off - ons;
        
        onsets{z} = ons;
        durations{z} = dur;
        names{z} = titleNamesShort{z,1};
        
    end
    
    %Insert instructions variable
    [tempOns,tempSort] = sort([onsets{1,1},onsets{1,2}]);
    onsets{3} = tempOns -3.3;
    durations{3} = ones(1,size(tempOns,2))*3.3;
    names{3} = 'Instructions';
    
    %Insert crosses variable
    onsets{4} = onsets{3} - 16;
    tempDur = [durations{1,1},durations{1,2}];
    tempDur2 = tempDur(tempSort);
    onsets{4}(1,end+1) = tempOns(1,end)+tempDur(1,end);
    durations{4} = ones(1,size(onsets{4},2))*16;
    names{4} = 'Crosses';
    
    %Remove start of scan
    for x = 1:size(onsets,2)
        onsets{x} = onsets{x} - onsets{4}(1,1);
    end
        
    save([pathName,ID,filesep,'Timing',filesep,ID,'_',cond{i,1},'.mat'],'onsets','durations','names');
end