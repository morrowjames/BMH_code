[a b] = uigetfile ('*.dat', 'MultiSelect','On');
%MRS_struct = GannetLoad({'meas_MID00078_FID43774_1_MEGA_GABA_HK.dat'},{'meas_MID00079_FID43775_1_MEGA_GABA_HK_Water.dat'});


for i = 1:2:length(a);

MRS_struct = GannetLoad({[b a{i}]},{[b a{i+1}]});
MRS_struct = GannetFit (MRS_struct);
[c d e] = fileparts (a{i});
mkdir([b d]);
%movefile ([b '*.mat'], [b d]);
f = dir([b 'MRSf*']);
movefile ([b f(end).name filesep 'MRS_struct.mat'], [b d]);
end




