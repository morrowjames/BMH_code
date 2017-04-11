function MRS_struct = GannetDiscernDatatype(filename,MRS_struct)
%Use the file ending to determine file type

lastchar=filename;

last2char=lastchar((end-1):end);
last4char=lastchar((end-3):end);

if(strcmpi(last2char,'.7'))
    MRS_struct.p.vendor = 'GE';
    MRS_struct.p.Reference_compound='H2O';
elseif(strcmpi(last4char,'SDAT'))
    MRS_struct.p.vendor = 'Philips';
    if(strcmp(last4char,'SDAT'))
       MRS_struct.p.spar_string='SPAR';
    else
        MRS_struct.p.spar_string='spar';
    end
elseif(strcmpi(last2char,'TA'))
    MRS_struct.p.vendor = 'Philips_data';
elseif(strcmpi(last2char,'DA'))
    MRS_struct.p.vendor = 'Siemens';
elseif(strcmpi(last4char,'.DAT'))
    MRS_struct.p.vendor = 'Siemens_twix';
else    
    error('Unrecognised filetype: should end .7 .SDAT or .RDA')
end
end
