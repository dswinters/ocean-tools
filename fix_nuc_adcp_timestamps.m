function adcp = fix_nuc_adcp_timestamps(adcp)

% Don't do anything on empty adcp structures
if isempty(adcp)
    return
end

% First pass - this should catch most
offset = abs(adcp.time-adcp.nuc_time); % offset between nuc/adcp clocks
offset_mean = nanmean(offset);
idx = find(offset >= 1.1*offset_mean); % timestamps offset by 10% more than average
if ~isempty(idx)
    idx_interp = setdiff(sort([idx-1,idx+1]),idx);
    [~,iu] = unique(idx_interp);
    adcp.nuc_time(idx) = interp1(idx_interp(iu),...
                                 adcp.nuc_time(idx_interp(iu)),...
                                 idx);
end

% Second pass - the first pass misses instances where timestamps are off twice in a row
idx = find(offset > 1.1*offset_mean);
if ~isempty(idx)
    idx_all = 1:length(adcp.time);
    for i = 1:length(idx)
        idx_pre = max(idx_all(offset < 1.1*offset_mean & idx_all < idx(i)));
        idx_post = min(idx_all(offset < 1.1*offset_mean & idx_all > idx(i)));
        adcp.nuc_time(idx(i)) = interp1([idx_pre, idx_post],...
                                        adcp.nuc_time([idx_pre, idx_post]),...
                                        idx(i));
    end
end
