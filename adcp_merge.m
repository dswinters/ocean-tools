%% adcp_merge_mat.m --- 
% Usage: A = adcp_merge(A,B)
% Description: Concat the adcp data in structure B with that in A
% Inputs: A: adcp data structure from rdradcp.m
%         B:  "
% Outputs: A: concatenated ADCP data structure
% Notes: This doesn't do any error checking. ADCP files to be merged must
%        have compatible numbers of bins/bin spacing etc.
% Author: Dylan Winters
% Created: September 12 2016

function A = adcp_merge(A,B)

A.mtime        = cat(2 , A.mtime        , B.mtime);
A.number       = cat(2 , A.number       , B.number);
A.pitch        = cat(2 , A.pitch        , B.pitch);
A.roll         = cat(2 , A.roll         , B.roll);
A.heading      = cat(2 , A.heading      , B.heading);
A.pitch_std    = cat(2 , A.pitch_std    , B.pitch_std);
A.roll_std     = cat(2 , A.roll_std     , B.roll_std);
A.heading_std  = cat(2 , A.heading_std  , B.heading_std);
A.depth        = cat(2 , A.depth        , B.depth);
A.temperature  = cat(2 , A.temperature  , B.temperature);
A.salinity     = cat(2 , A.salinity     , B.salinity);
A.pressure     = cat(2 , A.pressure     , B.pressure);
A.pressure_std = cat(2 , A.pressure_std , B.pressure_std);
A.east_vel     = cat(2 , A.east_vel     , B.east_vel);
A.north_vel    = cat(2 , A.north_vel    , B.north_vel);
A.vert_vel     = cat(2 , A.vert_vel     , B.vert_vel);
A.error_vel    = cat(2 , A.error_vel    , B.error_vel);
A.corr         = cat(3 , A.corr         , B.corr);
A.status       = cat(3 , A.status       , B.status);
A.intens       = cat(3 , A.intens       , B.intens);
A.bt_range     = cat(2 , A.bt_range     , B.bt_range);
A.bt_vel       = cat(2 , A.bt_vel       , B.bt_vel);
A.bt_corr      = cat(2 , A.bt_corr      , B.bt_corr);
A.bt_ampl      = cat(2 , A.bt_ampl      , B.bt_ampl);
A.bt_perc_good = cat(2 , A.bt_perc_good , B.bt_perc_good);
A.perc_good    = cat(3 , A.perc_good    , B.perc_good);

