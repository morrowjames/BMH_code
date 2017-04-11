clear

datapath='/Users/joshua_hendrikse/Documents/My_Documents/PhD/Ex_rTMS_study/Data/Active/'; %Experiment directory with subject folders in it
cd (datapath);

filter_str=''; fileselect; pause(0.1); %This will let you select subjects to run analysis on

%handknob_summary = [];
for bigK=1:length(FILE),

    hippocampus = [];
    Parietal = [];
    sma = [];
    allsummary = [];
    
    %fi=figure(bigK); set(fi,'Name',FILE{bigK})
    OutName = FILE{bigK};
    if isdir([datapath,'/',FILE{bigK}]),
        gopre = fullfile(datapath,'/',FILE{bigK},'/mri/MRS/pre'); %Sdirectory
        %cd (gohere);
        hippo_pre = dir([gopre '/meas_*hippocampus']);
        parietal_pre = dir([gopre '/meas_*Parietal']);
        sma_pre = dir([gopre '/meas_*SMA']);
        
        hippo_pre = fullfile(gopre, cell({hippo_pre.name}));
        parietal_pre = fullfile(gopre, cell({parietal_pre.name}));
        sma_pre = fullfile(gopre, cell({sma_pre.name}));
    end
    
    if isdir([datapath,'/',FILE{bigK}]),
        gopost = fullfile(datapath,'/',FILE{bigK},'/mri/MRS/post'); %Sdirectory
        %cd (gohere);
        hippo_post = dir([gopost '/meas_*hippocampus']);
        parietal_post = dir([gopost '/meas_*Parietal']);
        sma_post = dir([gopost '/meas_*SMA']);
        
        hippo_post = fullfile(gopost, cell({hippo_post.name}));
        parietal_post = fullfile(gopost, cell({parietal_post.name}));
        sma_post = fullfile(gopost, cell({sma_post.name}));
    end
    
  hippo_all = [hippo_pre, hippo_post];
  parietal_all = [parietal_pre, parietal_post];
  sma_all = [sma_pre, sma_post];
  
  all_all = [hippo_pre, hippo_post, parietal_pre, parietal_post, sma_pre, sma_post];
  
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