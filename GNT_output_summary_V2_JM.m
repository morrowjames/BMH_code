clear
datapath='/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/RAW_DATA/'; %Experiment directory with subject folders in it
cd (datapath);
addpath= '/gpfs/M2Home/projects/Monash076/Morrowj/MATLAB/scripts/';

filter_str=''; fileselect; pause(0.1); %This will let you select subjects to run analysis on

%handknob_summary = [];
for bigK=1:length(FILE),

    DLPFC = [];
    VIS = [];
    allsummary = [];
    
    %fi=figure(bigK); set(fi,'Name',FILE{bigK})
    OutName = FILE{bigK};
    if isdir([datapath,'/',FILE{bigK}]),
        gopre = fullfile(datapath,'/',FILE{bigK},'/MRI/MRS/PRE'); %Sdirectory
        %cd (gohere);
        DLPFC_pre = dir([gopre '/meas_*DLPFC']);
        Visual_pre = dir([gopre '/meas_*Visual']);
        
        DLPFC_pre = fullfile(gopre, cell({DLPFC_pre.name}));
        Visual_pre = fullfile(gopre, cell({Visual_pre.name}));
    end
    
    if isdir([datapath,'/',FILE{bigK}]),
        gopost = fullfile(datapath,'/',FILE{bigK},'/MRI/MRS/POST'); %Sdirectory
        %cd (gohere);
        DLPFC_post = dir([gopost '/meas_*DLPFC']);
        Visual_post = dir([gopost '/meas_*Visual']);
        
        DLPFC_post = fullfile(gopost, cell({DLPFC_post.name}));
        Visual_post = fullfile(gopost, cell({Visual_post.name}));
    end
    
  DLPFC_all = [DLPFC_pre, DLPFC_post];
  Visual_all = [Visual_pre, Visual_post];
  
  all_all = [DLPFC_pre, DLPFC_post, Visual_pre, Visual_post];
  
  for d = 1:length(all_all),
  load ([all_all{d}, '/MRS_struct.mat']);
  if isfield(MRS_struct.out, 'WaterArea')
      MRS_struct.out.WaterArea = MRS_struct.out.WaterArea;
      MRS_struct.out.GABAconciu = MRS_struct.out.GABAconciu;
      MRS_struct.out.WaterFitError = MRS_struct.out.WaterFitError;
      MRS_struct.out.GABAIU_Error_w = MRS_struct.out.GABAIU_Error_w;
  else
      MRS_struct.out.WaterArea = [99999];
      MRS_struct.out.GABAconciu = [99999];
      MRS_struct.out.WaterFitError = [99999];
      MRS_struct.out.GABAIU_Error_w = [99999];
  end
  
  %yogabagaba = [MRS_struct.out.GABAArea MRS_struct.out.CrArea MRS_struct.out.ChoArea MRS_struct.out.WaterArea MRS_struct.out.GABAconcCr MRS_struct.out.GABAconcCho MRS_struct.out.GABAconciu MRS_struct.out.GABAFitError MRS_struct.out.GABAFWHM MRS_struct.out.CrFWHMHz MRS_struct.out.FreqStdevHz MRS_struct.out.WaterFitError MRS_struct.out.WaterSNR MRS_struct.out.GABAsnr MRS_struct.out.GABAIU_Error_w MRS_struct.out.CrFitError MRS_struct.out.GABAIU_Error_cr];
  yogabagaba = [MRS_struct.out.GABAArea MRS_struct.out.CrArea MRS_struct.out.ChoArea MRS_struct.out.WaterArea MRS_struct.out.GABAconcCr MRS_struct.out.GABAconcCho MRS_struct.out.GABAconciu MRS_struct.out.GABAFitError MRS_struct.out.GABAFWHM MRS_struct.out.CrFWHMHz MRS_struct.out.FreqStdevHz MRS_struct.out.WaterFitError MRS_struct.out.GABAIU_Error_w MRS_struct.out.CrFitError MRS_struct.out.GABAIU_Error_cr];
  yogabagaba = yogabagaba';
  allsummary = ([allsummary, yogabagaba]);
  
  save('009_NR_all_summary')
  end
  
%   for d = 1:length(hippo_all),
%   load ([hippo_all{d}, '/MRS_struct.mat']);
%   if exist('MRS_struct.out.WaterArea', 'var')
%       MRS_struct.out.WaterArea = MRS_struct.out.WaterArea;
%       MRS_struct.out.GABAconciu = MRS_struct.out.GABAconciu;
%       MRS_struct.out.WaterFitError = MRS_struct.out.WaterFitError;
%       MRS_struct.out.GABAIU_Error_w = MRS_struct.out.GABAIU_Error_w;
%   else
%       MRS_struct.out.WaterArea = [99999];
%       MRS_struct.out.GABAconciu = [99999];
%       MRS_struct.out.WaterFitError = [99999];
%       MRS_struct.out.GABAIU_Error_w = [99999];
%   end
%   
%   %yogabagaba = [MRS_struct.out.GABAArea MRS_struct.out.CrArea MRS_struct.out.ChoArea MRS_struct.out.WaterArea MRS_struct.out.GABAconcCr MRS_struct.out.GABAconcCho MRS_struct.out.GABAconciu MRS_struct.out.GABAFitError MRS_struct.out.GABAFWHM MRS_struct.out.CrFWHMHz MRS_struct.out.FreqStdevHz MRS_struct.out.WaterFitError MRS_struct.out.WaterSNR MRS_struct.out.GABAsnr MRS_struct.out.GABAIU_Error_w MRS_struct.out.CrFitError MRS_struct.out.GABAIU_Error_cr];
%   yogabagaba = [MRS_struct.out.GABAArea MRS_struct.out.CrArea MRS_struct.out.ChoArea MRS_struct.out.WaterArea MRS_struct.out.GABAconcCr MRS_struct.out.GABAconcCho MRS_struct.out.GABAconciu MRS_struct.out.GABAFitError MRS_struct.out.GABAFWHM MRS_struct.out.CrFWHMHz MRS_struct.out.FreqStdevHz MRS_struct.out.WaterFitError MRS_struct.out.GABAIU_Error_w MRS_struct.out.CrFitError MRS_struct.out.GABAIU_Error_cr];
%   yogabagaba = yogabagaba';
%   hippocampus = ([hippocampus, yogabagaba]);
%   end
%   
%   for d = 1:length(parietal_all),
%   load ([parietal_all{d}, '/MRS_struct.mat']);
%   if exist('MRS_struct.out.WaterArea', 'var')
%       MRS_struct.out.WaterArea = MRS_struct.out.WaterArea;
%       MRS_struct.out.GABAconciu = MRS_struct.out.GABAconciu;
%       MRS_struct.out.WaterFitError = MRS_struct.out.WaterFitError;
%       MRS_struct.out.GABAIU_Error_w = MRS_struct.out.GABAIU_Error_w;
%   else
%       MRS_struct.out.WaterArea = [99999];
%       MRS_struct.out.GABAconciu = [99999];
%       MRS_struct.out.WaterFitError = [99999];
%       MRS_struct.out.GABAIU_Error_w = [99999];
%   %yogabagaba = [MRS_struct.out.GABAArea MRS_struct.out.CrArea MRS_struct.out.ChoArea MRS_struct.out.WaterArea MRS_struct.out.GABAconcCr MRS_struct.out.GABAconcCho MRS_struct.out.GABAconciu MRS_struct.out.GABAFitError MRS_struct.out.GABAFWHM MRS_struct.out.CrFWHMHz MRS_struct.out.FreqStdevHz MRS_struct.out.WaterFitError MRS_struct.out.WaterSNR MRS_struct.out.GABAsnr MRS_struct.out.GABAIU_Error_w MRS_struct.out.CrFitError MRS_struct.out.GABAIU_Error_cr];
%   yogabagaba = [MRS_struct.out.GABAArea MRS_struct.out.CrArea MRS_struct.out.ChoArea MRS_struct.out.WaterArea MRS_struct.out.GABAconcCr MRS_struct.out.GABAconcCho MRS_struct.out.GABAconciu MRS_struct.out.GABAFitError MRS_struct.out.GABAFWHM MRS_struct.out.CrFWHMHz MRS_struct.out.FreqStdevHz MRS_struct.out.WaterFitError MRS_struct.out.GABAIU_Error_w MRS_struct.out.CrFitError MRS_struct.out.GABAIU_Error_cr];
%   yogabagaba = yogabagaba';
%   Parietal = ([Parietal, yogabagaba]);
%   end

%%   If you really want to, you could create a structure that saves the output for all subjects
  summary(bigK).participant = OutName;
  summary(bigK).allsummary = allsummary;
%   summary(bigK).hippocampus = hippocampus;
%   summary(bigK).Parietal = Parietal;
  
end







%b = [MRS_struct.out.GABAArea MRS_struct.out.CrArea MRS_struct.out.ChoArea MRS_struct.out.WaterArea MRS_struct.out.GABAconcCr MRS_struct.out.GABAconcCho MRS_struct.out.GABAconciu MRS_struct.out.GABAFitError MRS_struct.out.GABAFWHM MRS_struct.out.GABAsnr MRS_struct.out.CrFWHMHz MRS_struct.out.FreqStdevHz MRS_struct.out.WaterFitError MRS_struct.out.WaterSNR MRS_struct.out.GABAIU_Error_w MRS_struct.out.CrFitError MRS_struct.out.GABAIU_Error_cr]