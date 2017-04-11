%mp_enable_stimulator turns the stimulator on or off. 
%It expects two input arguments:
%1: the com port of the stimulator you want
%2: numeric input of either a 0 (to disable), or a 1 (to enable).

function msg = mp_enable_stimulator(port,onoff)

if nargin<1 %if function called without arguments, assume switch on
    onoff=1;
    display('No input argument given. Assuming you want to switch the stimulator on');
end

if  ~(onoff==1 || onoff==0)
    display('you fail. type 1 to enable stimulator, 0 to disable it');
    return
end

front=['FE';'02'];
back=['FF'];

middle=['02'; '0' dec2hex(onoff)];

cs=mp_checksum_calculator(middle);
if length(cs)==1
    cs=['0' cs];
end

msg=hex2dec([front;middle;cs;back]);

fwrite(port,msg);

end