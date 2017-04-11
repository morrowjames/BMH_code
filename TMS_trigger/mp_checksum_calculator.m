%mp_checksum_calculator calculates CRC checksums for arbitrary polynomials
%first input argument should be a column of hex pairs such as ['1A';'2B';'3C'].
%second input argument should be the polynomial in binary form.
%if no second input argument is given the Dallas/Maxim CRC8 polynomial is
%used ie [1 0 0 1 1 0 0 0 1];

function sum=mp_checksum_calculator(in,poly)

if nargin<2 %ie, not polynomial specified
    poly=[1 0 0 1 1 0 0 0 1];
end

%take the hex bytes, convert them into binary (LSB first) and concatenate
msg=[];
for i=1:size(in,1)
    msg=[msg fliplr(str2num(dec2bin(hex2dec(in(i,:)),8)')')];
end
% the 'fliplr' is there because the system wants the least significant bit first
                                       
msg=[msg 0 0 0 0 0 0 0 0]; %append zeroes

count=length(poly);

reg=msg(1:count); %initialise

for b=(count+1):length(msg)
   if reg(1)==1;
   reg=xor(reg,poly); %xor msg w/ poly
   end
   reg=[reg(2:end) msg(b)]; %shift one place along and bring next digit down
end;

%NEEDS TO BE FLIPPED ONCE MORE!
if reg(1)==1;
reg=xor(reg,poly); %xor msg w/ poly
end

reg=fliplr(reg);
% again, least significant bit first

sum=dec2hex(bin2dec(num2str(reg(1:8))));
end

