%% adcp_merge_mat.m --- 
% Usage: A = adcp_merge_mat(A,files)
% Description: Load adcp data in FILES and concat it with adcp data in A
% Inputs: A: adcp data structure from rdradcp.m
%         files: cell array of strings w/ full paths to ADCP data files
% Outputs: A: concatenated ADCP data structure
% Notes: This doesn't do any error checking. ADCP files to be merged must
%        have compatible numbers of bins/bin spacing etc.
% Author: Dylan Winters
% Created: September 12 2016

function A = adcp_merge_mat(A,files)

adcp_varname = 'A'; % name of ADCP data structure in .mat files
i0 = 1;
if isempty(A)
    tmp = load(files{1},adcp_varname);
    A = tmp.(adcp_varname);
    i0 = 2;
end
for i = i0:length(files)
    disp(sprintf('%d entries, %d of %d files',length(A.mtime),i,length(files)))
    tmp = load(files{i},adcp_varname);
    A.mtime        = cat(2 , A.mtime        , tmp.(adcp_varname).mtime);
    A.number       = cat(2 , A.number       , tmp.(adcp_varname).number);
    A.pitch        = cat(2 , A.pitch        , tmp.(adcp_varname).pitch);
    A.roll         = cat(2 , A.roll         , tmp.(adcp_varname).roll);
    A.heading      = cat(2 , A.heading      , tmp.(adcp_varname).heading);
    A.pitch_std    = cat(2 , A.pitch_std    , tmp.(adcp_varname).pitch_std);
    A.roll_std     = cat(2 , A.roll_std     , tmp.(adcp_varname).roll_std);
    A.heading_std  = cat(2 , A.heading_std  , tmp.(adcp_varname).heading_std);
    A.depth        = cat(2 , A.depth        , tmp.(adcp_varname).depth);
    A.temperature  = cat(2 , A.temperature  , tmp.(adcp_varname).temperature);
    A.salinity     = cat(2 , A.salinity     , tmp.(adcp_varname).salinity);
    A.pressure     = cat(2 , A.pressure     , tmp.(adcp_varname).pressure);
    A.pressure_std = cat(2 , A.pressure_std , tmp.(adcp_varname).pressure_std);
    A.east_vel     = cat(2 , A.east_vel     , tmp.(adcp_varname).east_vel);
    A.north_vel    = cat(2 , A.north_vel    , tmp.(adcp_varname).north_vel);
    A.vert_vel     = cat(2 , A.vert_vel     , tmp.(adcp_varname).vert_vel);
    A.error_vel    = cat(2 , A.error_vel    , tmp.(adcp_varname).error_vel);
    A.corr         = cat(3 , A.corr         , tmp.(adcp_varname).corr);
    A.status       = cat(3 , A.status       , tmp.(adcp_varname).status);
    A.intens       = cat(3 , A.intens       , tmp.(adcp_varname).intens);
    A.bt_range     = cat(2 , A.bt_range     , tmp.(adcp_varname).bt_range);
    A.bt_vel       = cat(2 , A.bt_vel       , tmp.(adcp_varname).bt_vel);
    A.bt_corr      = cat(2 , A.bt_corr      , tmp.(adcp_varname).bt_corr);
    A.bt_ampl      = cat(2 , A.bt_ampl      , tmp.(adcp_varname).bt_ampl);
    A.bt_perc_good = cat(2 , A.bt_perc_good , tmp.(adcp_varname).bt_perc_good);
    A.perc_good    = cat(3 , A.perc_good    , tmp.(adcp_varname).perc_good);
end
