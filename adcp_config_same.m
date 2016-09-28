%%% adcp_config_same.m --- 
%% Usage: adcp_config_same(c1,c2)
%% Description: Check whether two ADCP configurations are the same
%% Inputs: c1,c2: 'config' fields from ADCP files read with rdradcp
%% Outputs: same: boolean, true if configs are the same
%% 
%% Author: Dylan Winters
%% Created: September 12 2016

function same = adcp_config_same(c1,c2,flds,varargin)
same = true;
verbose = false;
if ismember('-v',varargin); verbose=true;end


% Contents of each field identical?
for i = 1:length(flds)
    if ~all(c1.(flds{i}) == c2.(flds{i}))
        same = false;
        if verbose
            disp(sprintf('%s field not the same!',flds{i}))
        end
        return
    end
end

end