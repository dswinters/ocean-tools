%% adcp_get_serial.m
% Usage: serial = adcp_get_serial(A)
% Description: get the real serial number from an ADCP
%              data structure created by rdradcp.m.
% Inputs: A - rdradcp data structure
% Outputs: serial - ADCP's serial number
% 
% Author: Dylan Winters
% Created: 2017-03-24

function serial = adcp_get_serial(A);

% The ADCP's serial number is stored in 4 bytes.
% rdradcp returns the decimal representation of these bytes
% as a vector in the "config.remus_serialnum" field. We must
% convert them back to binary, concatenat them, and then convert
% them to decimal.
s = A.config.remus_serialnum;
b = [dec2bin(s(4),8) dec2bin(s(3),8) dec2bin(s(2),1) dec2bin(s(1),1)];
serial = bin2dec(b);
