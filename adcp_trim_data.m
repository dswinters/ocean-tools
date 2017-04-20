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
    [adcp(:).info] = deal({});
end

for i = 1:length(trim_methods)
    mes = 'Data removed ';
    switch trim_methods(i).name
      %% Do nothing
      case 'none'
        return

      %% Remove data below a percentage of bottom-track depth
      case 'BT'
        % Make sure the BT parameter is given as a percentage
        pct = trim_methods(i).params;
        if pct<1
            warning(sprintf(['BT percentage parameter was given as %.2f. ',...
                             'Assuming this means %.2f%%'],pct,100*pct));
            pct=100*pct;
        end
        mes = [mes sprintf('below %.2f%% of bottom-track depth',pct)];
        for ib = 1:adcp.config.n_beams;
            bd = adcp.bt_range(ib,:);
            [bd2 d2] = meshgrid(bd,adcp.config.ranges);
            bd(bd<5) = inf; % don't trust very shallow BT depths
            mask = d2>pct/100*bd2;
            vb = squeeze(adcp.vel(:,ib,:));
            vb(mask) = NaN;
            adcp.vel(:,ib,:) = vb;
        end
        adcp.info = cat(1,adcp.info,mes);

      %% Remove data using echo intensity edge-detection
      case 'ei_edge'
        mes = [mes 'using beam-specific echo intensity edge detection'];
        for ib = 1:adcp.config.n_beams
            edges = edge(squeeze(adcp.intens(:,ib,:)),'Sobel','horizontal');
            mask = (cumsum(edges) > 0);
            vb = squeeze(adcp.vel(:,ib,:));
            vb(mask) = NaN;
            adcp.vel(:,ib,:) = vb;
        end
        adcp.info = cat(1,adcp.info,mes);

      %% Remove data using correlation edge-detection
      case 'corr_edge'
        switch trim_methods(i).params
          case 'beam'
            mes = [mes 'using beam-specific correlation edge detection'];
            for ib = 1:adcp.config.n_beams
                edges = edge(squeeze(adcp.corr(:,ib,:)),'Sobel','horizontal');
                mask = (cumsum(edges) > 0);
                vb = squeeze(adcp.vel(:,ib,:));
                vb(mask) = NaN;
                adcp.vel(:,ib,:) = vb;
            end
          case 'avg'
            mes = [mes 'using beam-averaged correlation edge detection'];
            edges = edge(squeeze(nanmean(adcp.corr,2)),'Sobel','horizontal');
            mask = (cumsum(edges) > 0);
            for ib = 1:adcp.config.n_beams
                vb = squeeze(adcp.vel(:,ib,:));
                vb(mask) = NaN;
                adcp.vel(:,ib,:) = vb;
            end
        end
        adcp.info = cat(1,adcp.info,mes);

      %% Remove data below cutoff depth
      case 'cutoff'
        cutoff = trim_methods(i).params;
        mask = false(adcp.config.n_beams,length(adcp.mtime));
        mask(adcp.config.ranges > cutoff,:) = true;
        for ib=1:adcp.config.n_beams
            vb = squeeze(adcp.vel(:,ib,:));
            vb(mask) = NaN;
            adcp.vel(:,ib,:) = vb;
        end
        mes = [mes sprintf('below %.2fm depth',cutoff)];
        adcp.info = cat(1,adcp.info,mes);
        end        
    end
end

