%%% adcp_index.m --- 
%% Usage: A = adcp_index(A,idx)
%% Description: returns an ADCP data structure with only 
%%%             the specified time indices
%% Inputs: A - ADCP structure
%% Outputs: A - ADCP structure with limited time indices
%% 
%% Author: Dylan Winters
%% Created: August 10 2016

function A = adcp_index(A,idx)

A.mtime        = A.mtime(:,idx);
A.number       = A.number(:,idx);
A.pitch        = A.pitch(:,idx);
A.roll         = A.roll(:,idx);
A.heading      = A.heading(:,idx);
A.pitch_std    = A.pitch_std(:,idx);
A.roll_std     = A.roll_std(:,idx);
A.heading_std  = A.heading_std(:,idx);
A.depth        = A.depth(:,idx);
A.temperature  = A.temperature(:,idx);
A.salinity     = A.salinity(:,idx);
A.pressure     = A.pressure(:,idx);
A.pressure_std = A.pressure_std(:,idx);
A.east_vel     = A.east_vel(:,idx);
A.north_vel    = A.north_vel(:,idx);
A.vert_vel     = A.vert_vel(:,idx);
A.error_vel    = A.error_vel(:,idx);
A.corr         = A.corr(:,:,idx);
A.status       = A.status(:,:,idx);
A.intens       = A.intens(:,:,idx);
A.bt_range     = A.bt_range(:,idx);
A.bt_vel       = A.bt_vel(:,idx);
A.bt_corr      = A.bt_corr(:,idx);
A.bt_ampl      = A.bt_ampl(:,idx);
A.bt_perc_good = A.bt_perc_good(:,idx);
A.perc_good    = A.perc_good(:,:,idx);