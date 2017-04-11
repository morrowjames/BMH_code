%mp_start_train triggers a pulse train that has been preprogrammed into the
%stimulator. It takes no input arguments.

function msg = mp_start_train(port,junk)

if nargin>1 %if function called with arguments
   display('WARNING: this function only wants one input argument (the serial port id)');
end

front=['FE';'01'];
back=['FF'];

middle=['04'];

cs=mp_checksum_calculator(middle);
if length(cs)==1
    cs=['0' cs];
end

msg=hex2dec([front;middle;cs;back]);

fwrite(port,msg);

end