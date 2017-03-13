%% adcp_trim_data.m
% Usage: adcp_trim_data(adcp,trim_methods)
% Description: Remove data from bad depth cells
% Inputs:      adcp - adcp data structure
%      trim_methods - A struct with fields 'name' and 'params'
%                     -To remove data below a percent of BT depth:
%                        method.name = 'BT', method.params = percent
%                     -To remove data below a depth cutoff:
%                        method.name = 'cutoff', method.params = maxdepth
%                     -To remove data using echo intensity edge detection
%                        method.name = 'EI', method.params = [];
% Outputs: adcp - modified adcp data structure
% 
% Author: Dylan Winters
% Created: 2017-03-13

function [adcp] = adcp_trim_data(adcp,trim_methods)

if ~isfield(adcp,'info');
    adcp.info = {};
end

bnames = {'east','north','vert','error'};
for i = 1:length(trim_methods)
    mes = 'Data removed ';
    switch trim_methods(i).name
      %% Remove data below a percentage of bottom-track depth
      case 'BT'
        pct = trim_methods(i).params;
        mes = [mes 'below %.2f%% of bottom-track depth'];
        for ib = 1:4;
            bd = adcp.bt_range(ib,:);
            bd(bd<5) = inf; % don't trust very shallow BT depths
            [bd2 d2] = meshgrid(bd,adcp.config.ranges);
            adcp.([bnames{ib} '_vel'])(d2>pct/100*bd2) = NaN;
        end
        adcp.info = cat(1,adcp.info,mes);

      %% Remove data using echo intensity edge-detection
      case 'EI'
        mes = [mes 'using echo intensity edge detection'];
        % look for sharp increases in echo intensity, remove data below
        for ib = 1:4
            edges = edge(squeeze(adcp.intens(:,ib,:)),'Sobel','horizontal');
            bt_mask = (cumsum(edges) > 0);
            adcp.([bnames{ib} '_vel'])(bt_mask) = NaN;
        end
        adcp.info = cat(1,adcp.info,mes);

      %% Remove data below cutoff
      case 'cutoff'
        cutoff = trim_methods(i).params;
        mask = false(size(adcp.east_vel));
        mask(adcp.config.ranges > cutoff,:) = true;
        for ib=1:4
            adcp.([bnames{ib} '_vel'])(mask) = NaN;
        end
        mes = [mes sprintf('below %.2fm depth',cutoff)];
        adcp.info = cat(1,adcp.info,mes);
        end        
    end
end

