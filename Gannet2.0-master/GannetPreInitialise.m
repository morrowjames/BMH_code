function MRS_struct=GannetPreInitialise(MRS_struct)

% Some of these parameters will be overwritten by correct values are stored
% in the data headers.


%Acquisition Parameters
    MRS_struct.p.sw=2000;  % sw taken from header for all formats except Philips .data
   % MRS_struct.p.sw=5000; %This should be parsed from headers where possible
    MRS_struct.p.npoints=2080; %This is twice the acquired points for TWIX data;
    %This should be parsed from headers where possible
    MRS_struct.p.TR=2000;%This should be parsed from headers where possible
    MRS_struct.p.TE=68; %This should be parsed from headers where possible
    MRS_struct.p.LarmorFreq=127; %This should be parsed from headers where possible
    %In general, LarmorFreq is 127.8 on Philips,
    MRS_struct.p.Nwateravg = 16; %Needed for GE
    MRS_struct.p.target='GABA'; %Other options are 'GSH' and 'Glx' (now also implemented) 'GABAGlx' in dev
    MRS_struct.p.ONOFForder='onfirst';
    %Options are MRS_struct.ONOFForder='onfirst' or 'offfirst';
    MRS_struct.p.Water_Positive=1; %For Philips MOIST ws, set to 0.
    %Siemens header information differs between versions
    %switch for different versions
    MRS_struct.p.Siemens_type = 1; %1 = TIM TRIO WIP 2 = Near seq 3 =Skyra WIP; 4=Prisma (VD13C); 5=Prisma(Minnesota); 6 = Jamie's VE11B (Jena); 7 = Minnesota (VE11C)
    
%Analysis Parameters
    MRS_struct.p.LB = 3;
    MRS_struct.p.ZeroFillTo = 32768;
    %AlignTo planned options: Cr; Cho; NAA; H20; CrOFF
    MRS_struct.p.AlignTo = 'SpecReg'; %SpecReg default and recommended
    
%Output Parameters
    MRS_struct.p.mat = 1; %1 = YES, save MRS_struct as .mat file.
    MRS_struct.p.sdat = 0; %1 = YES, save MRS_struct as .sdat file.
end
