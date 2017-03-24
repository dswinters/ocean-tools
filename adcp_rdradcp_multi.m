%% adcp_rdradcp_multi.m --- 
% Usage: adcp_rdradcp_multi(fn)
% Description: loads and concatenates data from multiple binary ADCP files
% Inputs: fn - cell array of file names
% Outputs: A - ACDP data structure
% 
% Author: Dylan Winters
% Created: August 10 2016

function A = adcp_rdradcp_multi(fn,skip)

A = [];
ng = 0;
for i = 1:length(fn)
    adcp = [];
    try
        adcp = rdradcp(fn{i},skip,-1);
        adcp.file = 0*adcp.mtime + i;
        ng = ng+1;
    catch err
    end

    % Concat data from all files
    if ng == 1;
        A = adcp;
    else
        if ~isempty(adcp)
            A.mtime        = cat(2 , A.mtime        , adcp.mtime);
            A.number       = cat(2 , A.number       , adcp.number);
            A.pitch        = cat(2 , A.pitch        , adcp.pitch);
            A.roll         = cat(2 , A.roll         , adcp.roll);
            A.heading      = cat(2 , A.heading      , adcp.heading);
            A.pitch_std    = cat(2 , A.pitch_std    , adcp.pitch_std);
            A.roll_std     = cat(2 , A.roll_std     , adcp.roll_std);
            A.heading_std  = cat(2 , A.heading_std  , adcp.heading_std);
            A.depth        = cat(2 , A.depth        , adcp.depth);
            A.temperature  = cat(2 , A.temperature  , adcp.temperature);
            A.salinity     = cat(2 , A.salinity     , adcp.salinity);
            A.pressure     = cat(2 , A.pressure     , adcp.pressure);
            A.pressure_std = cat(2 , A.pressure_std , adcp.pressure_std);
            A.east_vel     = cat(2 , A.east_vel     , adcp.east_vel);
            A.north_vel    = cat(2 , A.north_vel    , adcp.north_vel);
            A.vert_vel     = cat(2 , A.vert_vel     , adcp.vert_vel);
            A.error_vel    = cat(2 , A.error_vel    , adcp.error_vel);
            A.corr         = cat(3 , A.corr         , adcp.corr);
            A.status       = cat(3 , A.status       , adcp.status);
            A.intens       = cat(3 , A.intens       , adcp.intens);
            A.bt_range     = cat(2 , A.bt_range     , adcp.bt_range);
            A.bt_vel       = cat(2 , A.bt_vel       , adcp.bt_vel);
            A.bt_corr      = cat(2 , A.bt_corr      , adcp.bt_corr);
            A.bt_ampl      = cat(2 , A.bt_ampl      , adcp.bt_ampl);
            A.bt_perc_good = cat(2 , A.bt_perc_good , adcp.bt_perc_good);
            A.perc_good    = cat(3 , A.perc_good    , adcp.perc_good);
            A.file         = cat(2 , A.file         , adcp.file);
        end
    end
end

if isempty(A)
    return
end

% Remove incomplete ensembles
nt = min([size(A.mtime,2);
          size(A.number,2);
          size(A.pitch,2);
          size(A.roll,2);
          size(A.heading,2);
          size(A.pitch_std,2);
          size(A.roll_std,2);
          size(A.heading_std,2);
          size(A.depth,2);
          size(A.temperature,2);
          size(A.salinity,2);
          size(A.pressure,2);
          size(A.pressure_std,2);
          size(A.east_vel,2);
          size(A.north_vel,2);
          size(A.vert_vel,2);
          size(A.error_vel,2);
          size(A.corr,3);
          size(A.status,3);
          size(A.intens,3);
          size(A.bt_range,2);
          size(A.bt_vel,2);
          size(A.bt_corr,2);
          size(A.bt_ampl,2);
          size(A.bt_perc_good,2);
          size(A.perc_good,3)]);
kp = A.mtime>0;
A.mtime        = A.mtime(:,kp(1:nt));
A.number       = A.number(:,kp(1:nt));
A.pitch        = A.pitch(:,kp(1:nt));
A.roll         = A.roll(:,kp(1:nt));
A.heading      = A.heading(:,kp(1:nt));
A.pitch_std    = A.pitch_std(:,kp(1:nt));
A.roll_std     = A.roll_std(:,kp(1:nt));
A.heading_std  = A.heading_std(:,kp(1:nt));
A.depth        = A.depth(:,kp(1:nt));
A.temperature  = A.temperature(:,kp(1:nt));
A.salinity     = A.salinity(:,kp(1:nt));
A.pressure     = A.pressure(:,kp(1:nt));
A.pressure_std = A.pressure_std(:,kp(1:nt));
A.east_vel     = A.east_vel(:,kp(1:nt));
A.north_vel    = A.north_vel(:,kp(1:nt));
A.vert_vel     = A.vert_vel(:,kp(1:nt));
A.error_vel    = A.error_vel(:,kp(1:nt));
A.corr         = A.corr(:,:,kp(1:nt));
A.status       = A.status(:,:,kp(1:nt));
A.intens       = A.intens(:,:,kp(1:nt));
A.bt_range     = A.bt_range(:,kp(1:nt));
A.bt_vel       = A.bt_vel(:,kp(1:nt));
A.bt_corr      = A.bt_corr(:,kp(1:nt));
A.bt_ampl      = A.bt_ampl(:,kp(1:nt));
A.bt_perc_good = A.bt_perc_good(:,kp(1:nt));
A.perc_good    = A.perc_good(:,:,kp(1:nt));
A.file         = A.file(:,kp);

% Include which files were read
A.files = cell(length(fn),1);
for i = 1:length(fn)
    [fdir fname fext] = fileparts(fn{i});
    A.files{i} = [fname fext];
end
