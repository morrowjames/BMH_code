%mp_single_trigger sends a single trigger pulse via the COM port to the
%stimulator. it takes one input argument; the port object handle.

function msg = mp_single_trigger(port,junk)

if nargin>1 %if function called with arguments
   display('WARNING: this function only wants one input argument (the serial port id)');
end

front=['FE';'02'];
back=['FF'];

middle=['03';'01'];

cs=mp_checksum_calculator(middle);
if length(cs)==1
    cs=['0' cs];
end

msg=hex2dec([front;middle;cs;back]);

fwrite(port,msg);

end