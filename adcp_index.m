%% adcp_index.m --- 
% Usage: A = adcp_index(A,idx)
% Description: returns an ADCP data structure with only 
%              the specified time indices. This function
%              assumes that the time dimension is last for 
%              all fields in the ADCP data structure.
% Inputs: A - ADCP structure
% Outputs: A - ADCP structure with limited time indices
% 
% Author: Dylan Winters
% Created: August 10 2016

function A = adcp_index(A,idx)

flds = setdiff(fields(A),{'config','files','info','masks','field_info','gps'});
for i = 1:length(flds)
    nd = ndims(A.(flds{i}));
    if nd == 2 & ~isstr(A.(flds{i}));
        A.(flds{i}) = A.(flds{i})(:,idx);
    elseif nd == 3 & ~isstr(A.(flds{i}));
        A.(flds{i}) = A.(flds{i})(:,:,idx);
    end
end
