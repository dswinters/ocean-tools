function adcp = adcp_filter(adcp,flt)

if ~isfield(flt,'name')
    return
end

if ~isfield(adcp,'info')
    [adcp.info(:)] = deal({});
end

switch flt.name
  case 'rotmax' % Remove data where rate of turn exceeds some value
    mes = sprintf('Data removed where rate of turn > %.2f deg/s',flt.params);
    dn = adcp.mtime(1:end-1) + diff(adcp.mtime)/2;
    dh = diff(adcp.heading);
    dh(dh>180) = dh(dh>180) - 360;
    dh(dh<-180) = dh(dh<-180) + 360;
    dt = diff(adcp.mtime)*86400;
    rot = interp1(dn,dh./dt,adcp.mtime);
    idx = rot > flt.params;
    adcp.vel(:,:,idx) = NaN;
    adcp.info = cat(1,adcp.info,mes);
  case 'corrmin'
    mes = sprintf('Data removed where correlation < %d counts',flt.params);
    adcp.vel(adcp.corr<flt.params) = NaN;
    adcp.info = cat(1,adcp.info,mes);
end







