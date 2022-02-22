%% nav_ltln2vel.m
% Usage: nav_ltln2vel(lt,ln,dn)
% Description: Convert the lat and lon measurements in LT and 
%              LN to an xy grid centered at the mean location.
%              Create a timeseries of x&y velocity based on the
%              timestamps in DN.
% Inputs: lt - latitude (degrees north)
%         ln - longitude (degrees east)
%         dn - matlab datenum
% Outputs: vx,vy - velocities, east & north
% 
% Author: Dylan Winters
% Created: 2016-09-16

function [vx, vy] = gps_ltln2vel(lt,ln,dn)

% remove non-unique timestamps
dn0 = dn;
[~,idx] = unique(dn);
dn = dn(idx);
lt = lt(idx);
ln = ln(idx);
% remove NaNs
idx = ~isnan(dn.*lt.*ln);
dn = dn(idx);
lt = lt(idx);
ln = ln(idx);

% convert to xy (need at least 2 non-nan points)
if sum(idx)>2
    wgs84 = referenceEllipsoid('wgs84','m');
    lt0 = nanmean(lt);
    ln0 = nanmean(ln);
    lt2y = distance('rh',lt0-0.5,ln0,lt0+0.5,ln0,wgs84);
    ln2x = distance('rh',lt0,ln0-0.5,lt0,ln0+0.5,wgs84);
    % lt2y = abs(40000000/360) ;   % meters N/S per degree N
    % ln2x = scl*cosd(lt0)     ;   % meters E/W per degree W at latitude lt0
    y  =  lt2y * (lt-lt0) ;    % meters N/S
    x  =  ln2x * (ln-ln0) ;    % meters E/W
    dt = diff(dn)*86400;
    t  = dn(1:end-1) + diff(dn)/2;
    vx = interp1(t, diff(x)./dt, dn0,'linear','extrap');
    vy = interp1(t, diff(y)./dt, dn0,'linear','extrap');
else
    vx = nan*dn0;
    vy = nan*dn0;
end
