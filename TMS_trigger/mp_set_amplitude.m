%mp_set_amplitude askes for a human-friendly amplitue value and generates a
%Rapid-friendly message to set the stimulator intensity to that value.
%first input argument should be the id of the (already opened) serial port.
%second input argument should be desired percentage intensity in decimal form.

function msg = mp_set_amplitude(port,amp)

if ~isnumeric(amp) || rem(amp,1)~=0 || amp<0 || amp>100
    display('you fail. set the amplitude as an integer of range 0-100');
    return
end

front=['FE';'02'];
back=['FF'];

d=dec2hex(amp);
if length(d)==1
    d=['0' d];
end

middle=['01';d];

cs=mp_checksum_calculator(middle);
if length(cs)==1
    cs=['0' cs];
end

msg=hex2dec([front;middle;cs;back]);

fwrite(port,msg);
end