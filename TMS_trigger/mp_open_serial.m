%mp_open_serial opens a named serial (COM) port. It takes one argument,
%which is a string naming this port (eg, 'COM2'). If called without
%arguments it will open COM1.

function s=mp_open_serial(com_port)

if nargin<1 %if function called without arguments
   com_port='COM1';
   display('no COM port specified. Opening COM1 (default)');
end

s = serial(com_port);
s.BaudRate = 38400;
s.DataBits = 8;
s.Parity = 'none';
s.StopBits = 1;
s.InputBufferSize= 8;
s.OutputBufferSize= 7;

fopen(s)

end