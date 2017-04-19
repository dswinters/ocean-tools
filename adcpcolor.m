%% adcpcolor.m
% Usage: adcpcolor(varname,beam)
% Description: Create a composite pcolor plot showing dual-profile
%              ADCP data in the same axes.
% Inputs: varname - variable to use for color values (string)
%            beam - beam number to use for plotting data with a beam
%                   dimension (e.g. velocity or echo intensity)
% Outputs: --
% 
% Author: Dylan Winters
% Created: 2017-04-17

function pcols = adcpcolor(adcp,varname,beam);


if length(adcp) > 1;
% Insert some slices of NaN's between every profile switch.
% This prevents the pcolor function from stretching data to fill
% in gaps, and instead creates blank spaces between sets of pings.
flds = {varname,'mtime'};
for i = 1:length(adcp)
    nb = adcp(i).config.n_beams;
    nc = adcp(i).config.n_cells;
    nt = length(adcp(i).mtime);
    dt = diff(adcp(i).mtime);
    prof_switched = find([false, dt>nanmean(dt)]);
    pidx = length(prof_switched);
    while pidx>1
        for f = 1:length(flds)
            % Make 2d fields temporarily 3d
            is2d = ndims(adcp(i).(flds{f})) == 2;
            if is2d
                adcp(i).(flds{f}) = permute(adcp(i).(flds{f}),[1 3 2]);
            end
            % Insert a slice of NaN's in the time dimension
            preidx = 1:prof_switched(pidx)-1;
            postidx = prof_switched(pidx):nt;
            nans = nan(size(adcp(i).(flds{f})(:,:,1)));
            if strcmp('mtime',flds{f});
                nans = adcp(i).mtime(preidx(end)) + ...
                       diff(adcp(i).mtime(preidx(end-1:end)));
            end
            adcp(i).(flds{f}) = cat(3,...
                                    adcp(i).(flds{f})(:,:,preidx([1:end])),...
                                    nans,...
                                    adcp(i).(flds{f})(:,:,postidx));
            if is2d
                adcp(i).(flds{f}) = squeeze(adcp(i).(flds{f}))';
            end
        end
        pidx = pidx-1;
        nt = length(adcp(i).mtime);
    end
end

%% Draw pcolor
for i = 1:length(adcp)
    cvals = squeeze(nanmean(adcp(i).(varname)(:,beam,:),2));
    pcols(i) = pcolor(adcp(i).mtime,...
                      adcp(i).config.ranges,...
                      cvals);
    hold on
    shading flat
end
hold off
else
    cvals = squeeze(nanmean(adcp.(varname)(:,beam,:),2));
    pcols = pcolor(adcp.mtime,...
                   adcp.config.ranges,...
                   cvals);
    shading flat
end