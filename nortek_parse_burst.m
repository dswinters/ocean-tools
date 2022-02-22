function burst = nortek_parse_burst(d)

getdat =@(dat,idx,type) double(typecast(reshape([dat(:,idx)]',[],1),type))';

burst = struct();
offset = getdat(d(1,:),2,'uint8'); % this is where vel/amp/etc start
config          = dec2bin(getdat(d,3:4,'uint16'),16);
burst.config.has_vel = logical(bin2dec(config(1,11)));
burst.config.has_amp = logical(bin2dec(config(1,10)));
burst.config.has_cor = logical(bin2dec(config(1,9)));
burst.config.has_alt = logical(bin2dec(config(1,8)));
burst.config.has_alt_raw = logical(bin2dec(config(1,7)));
burst.config.has_ast = logical(bin2dec(config(1,6)));
burst.config.has_ecs = logical(bin2dec(config(1,5)));
burst.config.has_ahrs = logical(bin2dec(config(1,4)));
burst.config.has_pg = logical(bin2dec(config(1,3)));
burst.config.has_std = logical(bin2dec(config(1,2)));
burst.valid.pressure = logical(bin2dec(config(:,16)));
burst.valid.temp = logical(bin2dec(config(:,15)));
burst.valid.compass = logical(bin2dec(config(:,14)));
burst.valid.tilt = logical(bin2dec(config(:,13)));
burst.serial    = getdat(d(1,:),5:8,'uint32');

% Clock
year            = getdat(d,9,'uint8') + 1900;
month           = getdat(d,10,'uint8') + 1;
day             = getdat(d,11,'uint8');
hour            = getdat(d,12,'uint8');
mins            = getdat(d,13,'uint8');
sec             = getdat(d,14,'uint8');
sec             = sec + getdat(d,15:16,'uint16')*1e-4;
burst.time      = datenum([year' month' day' hour' mins' sec'])';
burst.sos       = getdat(d,17:18,'uint16')/10; % m/s
burst.temp      = getdat(d,19:20,'int16')/100; % degrees
burst.pressure  = getdat(d,21:24,'int32')/1000; % dbar
burst.heading   = getdat(d,25:26,'uint16')/100; % Deg
burst.pitch     = getdat(d,27:28,'int16')/100; % Deg
burst.roll      = getdat(d,29:30,'int16')/100; % Deg
beaminfo        = dec2bin(getdat(d,31:32,'uint16'),16);
burst.ncells    = bin2dec(beaminfo(1,end-9:end));
switch beaminfo(1,end-11:end-10)
  case '00'; burst.coords = 'ENU';
  case '01'; burst.coords = 'XYZ';
  case '10'; burst.coords = 'BEAM';
end
burst.nbeams    = bin2dec(beaminfo(1,end-15:end-12));
burst.cell_size = getdat(d,33:34,'uint16')/1000; % m
burst.blanking  =  getdat(d,35:36,'uint16')/100; % m
burst.cor_nom  =  getdat(d,37,'uint8'); % percent
burst.pres_temp =  getdat(d,38,'uint8')/5-4; % deg c
burst.bat_voltage  = getdat(d,39:40,'uint16')/10; % volts
magx = getdat(d,41:42,'int16');
magy = getdat(d,43:44,'int16');
magz = getdat(d,45:46,'int16');
accx = getdat(d,47:48,'int16');
accy = getdat(d,49:50,'int16');
accz = getdat(d,51:52,'int16');
burst.raw_mag = [magx', magy', magz'];
burst.raw_acc = [accx', accy', accz'];

% TODO missing some fields around here

nb = burst.nbeams;
nc = burst.ncells;
nt = length(burst.time);
nd = nb*nc;

% Beams aren't necessarily in order
bnums = dec2bin(getdat(d(1,:),55:56,'uint16'),16);
bnums = [bin2dec(bnums(13:16)),...
         bin2dec(bnums(9:12)),...
         bin2dec(bnums(5:8)),...
         bin2dec(bnums(1:4))];
[~,bidx] = sort(bnums); % permutation to put beams in order
resh = @(x) permute(reshape(x,[nc nb nt]), [1 2 3]);


% Extract velocity, amplitude, and correlation
vel_scaling = getdat(d(1,:),59,'int8');
if burst.config.has_vel
    vel = getdat(d,offset + [1:2*nd],'int16') * 10^vel_scaling(1);
    burst.vel = resh(vel);
    burst.vel = burst.vel(:,bidx,:);
    offset = offset + 2*nd;
end
if burst.config.has_amp
    amp = getdat(d,offset + [1:nd],  'uint8');
    burst.amp = resh(amp);
    burst.amp = burst.amp(:,bidx,:);
    offset = offset + nd;
end
if burst.config.has_cor
    cor = getdat(d,offset + [1:nd],  'uint8');
    burst.cor = resh(cor);
    burst.cor = burst.cor(:,bidx,:);
    offset = offset + nd;
end

% More fields that may or may not exisxt depending on the config
if burst.config.has_alt
    burst.alt_distance = getdat(d,offset + [1:4],'single');
    burst.alt_quality = getdat(d,offset + [5:6],'uint16');
    % burst.alt_status = dec2bin(getdat(d,offset+[7:8],'uint16'));
    offset = offset + 8;
end

if burst.config.has_ast
    burst.ast_distance = getdat(d, offset + [1:4],'single');
    burst.ast_quality = getdat(d,offset + [5:6],'uint16');
    burst.ast_offset100us = getdat(d,offset + [7:8],'int16');
    burst.ast_pressure = getdat(d,offset + [9:12],'single');
    offset = offset + 20;
end

if burst.config.has_alt_raw
    % TODO
    offset = offset + 8;
end

if burst.config.has_ecs
    burst.ecs = 0.01*reshape(getdat(d,offset+[1:2*nc],'uint16'),nc,[]);
    offset = offset + 2*nc;
end

if burst.config.has_ahrs
    rotm = getdat(d, offset + [1:4*9],'single');
    offset = offset + 5*9;
    burst.ahrs_rotm = pagetranspose(reshape(rotm,3,3,[]));
    burst.ahrs_gyro = reshape(getdat(d, offset+[1:3*4],'single'),3,[])';
end

if burst.config.has_pg
    % TODO
    offset = offset + nc;
end

if burst.config.has_std
    % TODO
end
