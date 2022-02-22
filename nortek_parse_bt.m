function bt = nortek_parse_bt(d)

    getdat =@(dat,idx,type) double(typecast(reshape([dat(:,idx)]',[],1),type))';

    bt = struct();
    offset = getdat(d(1,:),2,'uint8');
    config = fliplr(dec2bin(getdat(d,3:4,'uint16'),16));
    bt.config.has_vel = logical(bin2dec(config(1,6)));
    bt.config.has_dist = logical(bin2dec(config(1,9)));
    bt.config.has_fom = logical(bin2dec(config(1,10)));
    bt.valid.pressure = logical(bin2dec(config(:,1)));
    bt.valid.temp = logical(bin2dec(config(:,2)));
    bt.valid.compass = logical(bin2dec(config(:,3)));
    bt.valid.tilt = logical(bin2dec(config(:,4)));

    % Clock
    year    = getdat(d,9,'uint8') + 1900;
    month   = getdat(d,10,'uint8') + 1;
    day     = getdat(d,11,'uint8');
    hour    = getdat(d,12,'uint8');
    mins    = getdat(d,13,'uint8');
    sec     = getdat(d,14,'uint8');
    sec     = sec + getdat(d,15:16,'uint16')*1e-4;
    bt.time = datenum([year' month' day' hour' mins' sec'])';
    bt.sos       = getdat(d,17:18,'uint16')/10; % m/s
    bt.temp      = getdat(d,19:20,'int16')/100; % degrees
    bt.pressure  = getdat(d,21:24,'int32')/1000; % dbar
    bt.heading   = getdat(d,25:26,'uint16')/100; % Deg
    bt.pitch     = getdat(d,27:28,'int16')/100; % Deg
    bt.roll      = getdat(d,29:30,'int16')/100; % Deg
    beaminfo     = dec2bin(getdat(d,31:32,'uint16'),16);

    bt.ncells    = bin2dec(beaminfo(1,end-9:end));
    switch beaminfo(1,end-11:end-10)
      case '00'; bt.coords = 'ENU';
      case '01'; bt.coords = 'XYZ';
      case '10'; bt.coords = 'BEAM';
    end
    bt.nbeams    = bin2dec(beaminfo(1,end-15:end-12));
    bt.cell_size = getdat(d,33:34,'uint16')/1000; % m
    bt.blanking  =  getdat(d,35:36,'uint16')/1000; % m
    bt.cor_nom  =  getdat(d,37,'uint8'); % percent

    % Beams aren't necessarily in order
    bnums = dec2bin(getdat(d(1,:),57:58,'uint16'),16);
    bnums = [bin2dec(bnums(13:16)),...
             bin2dec(bnums(9:12)),...
             bin2dec(bnums(5:8)),...
             bin2dec(bnums(1:4))];
    [~,bidx] = sort(bnums); % permutation to put beams in order
    resh = @(x) permute(reshape(x,[nc nb nt]), [1 2 3]);

    vel_scale = getdat(d(1,:),61,'int8');
    nb = bt.nbeams;
    if bt.config.has_vel
        vel = 10^vel_scale * getdat(d,offset + [1:4*nb],'int32');
        bt.vel = reshape(vel,nb,[]);
        offset = offset + 4*nb;
    end

    if bt.config.has_dist
        dist = getdat(d,offset + [1:4*nb],'int32')/1e6; % m
        bt.dist = 1e3*reshape(dist,nb,[]);
        offset = offset + 4*nb;
    end

    if bt.config.has_fom
        fom = getdat(d,offset + [1:2*nb],'uint16');
        bt.fom = 1e3*reshape(dist,nb,[]);
        offset = offset + 2*nb;
    end
