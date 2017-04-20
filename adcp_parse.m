function adcp = adcp_parse(files,varargin)

%% Setup
setup = struct();
if isstr(files)
    files = {files};
end
setup.files = files;
setup.nprofs = 1;

% ROSS computers insert an extra 7 bytes of data
% after each ensemble start indicator.
setup.ross_bytes = 0;
if ismember('ross',lower(varargin))
    setup.ross_bytes = 7;
end

dat = load_data(files);               % load all data
[h nbytes] = find_headers(dat,setup); % find header locations
setup.nens = length(h);
adcp = [];                            % initialize ADCP structure
aidx = [];


%% Fill adcp data structure with ensemble data
hidx = 1;
disp('Parsing')
while hidx < length(h);
    ensidx = h(hidx) + [0:nbytes(hidx)-1]; % Compute indices of ensemble data
    ensdata = double(dat(ensidx));              % Extract ensemble data
    [adcp aidx] = get_data(adcp,ensdata,aidx,setup); % Fill adcp structure with data
    hidx = hidx + 1;
    print_progress(hidx,length(h),0);
end


%% Remove empty ensembles
for i = 1:length(adcp)
    adcp(i) = rm_ens(adcp(i),isnan(adcp(i).mtime));
end


%%% Sub-functions

%----------------------------------------------------------
% Read and concatenate data in the given list of files
%----------------------------------------------------------
function dat = load_data(files)
dat = [];
for i = 1:length(files)
    if ~exist(files{i},'file')
        error('File not found: %s',files{i})
    end
    [~,fname,fext] = fileparts(files{i});
    disp(['Opening ' fname fext]);
    fd = fopen(files{i},'r','ieee-le');
    dat = cat(1,dat,uint8(fread(fd,inf,'uint8')));
    fclose(fd);
end

%----------------------------------------------------------
% Transform a sequence of n-bit integers into the single
% decimal number represented by their concatenated binary 
% sequence. For example, intcat([127 127],8) gives the 
% result of converting the sequence to 8-bit binary
% representation:
%   01111111 01111111
% concatenating it:
%   0111111101111111
% and converting it back to decimal:
%  >> 32639
% 
% Each row of the input matrix represents a separate 
% sequence.
%----------------------------------------------------------
function n = intcat(intseq,nbits)
n = double(intseq) * [2.^(nbits*[0:size(intseq,2)-1])]';

%----------------------------------------------------------
% Convert an unsigned integer to a signed integer
%----------------------------------------------------------
function n = uint2int(n,base)
m = 2^base/2;
n(n>m) = -2*m + n(n>m);

%----------------------------------------------------------
% Find all header start indices in raw data
%----------------------------------------------------------
function [h nbytes] = find_headers(dat,setup)

rb = setup.ross_bytes;
h = find(dat(1:end-1)==127 & dat(2:end)==127);

if length(h)==0
    error('No headers found!')
end

if h(end)+2 > length(dat)
    h = h(1:end-1);
end
nbytes = intcat([dat(h+2+rb),...
                 dat(h+3+rb)],8) + rb;

% Remove headers that appear mid-ensemble
i = 1;
while i<length(h)
    if h(i+1) < h(i)+nbytes(i)+2;
        h(i+1) = [];
        nbytes(i+1) = [];
    else
        i = i + 1;
    end
end

% Remove incomplete ensemble at end
if h(end) + nbytes(end) > length(dat)
    h = h(1:end-1);
    nbytes = nbytes(1:end-1);
end

% Remove ensembles with anomalous byte counts
[counts,unb] = hist(nbytes,unique(nbytes));
for i = 1:length(unb)
    if counts(i)/length(nbytes) < 0.1;
        h = h(nbytes~=unb(i));
        nbytes = nbytes(nbytes~=unb(i));
    end
end

%----------------------------------------------------------
% Initialize adcp data structure
%----------------------------------------------------------
function adcp = adcp_init(ensdata,setup)

adcp = struct();
adcp.config = get_fixed_leader(ensdata,setup);

nt = floor(setup.nens);
nc = adcp.config.n_cells;
nb = adcp.config.n_beams;
%
adcp.mtime        = nan(1,nt);
if setup.ross_bytes > 0
    adcp.ross_mtime = nan(1,nt);
end
adcp.number       = nan(1,nt);
adcp.pitch        = nan(1,nt);
adcp.roll         = nan(1,nt);
adcp.heading      = nan(1,nt);
adcp.pitch_std    = nan(1,nt);
adcp.roll_std     = nan(1,nt);
adcp.heading_std  = nan(1,nt);
adcp.depth        = nan(1,nt);
adcp.salinity     = nan(1,nt);
adcp.temperature  = nan(1,nt);
adcp.pressure     = nan(1,nt);
adcp.pressure_std = nan(1,nt);
adcp.vel          = nan(nc,nb,nt);
adcp.corr         = nan(nc,nb,nt);
adcp.status       = nan(nc,nb,nt);
adcp.intens       = nan(nc,nb,nt);
adcp.perc_good    = nan(nc,nb,nt);
adcp.bt_range     = nan(nb,nt);
adcp.bt_vel       = nan(nb,nt);
adcp.bt_corr      = nan(nb,nt);
adcp.bt_ampl      = nan(nb,nt);
adcp.bt_perc_good = nan(nb,nt);
adcp.files        = cell(length(setup.files),1);

for nf = 1:length(setup.files)
    [~,fname,fext] = fileparts(setup.files{nf});
    adcp.files{nf} = [fname fext];
end


%----------------------------------------------------------
% Return header information from the given ensemble.
%----------------------------------------------------------
function header = get_header(ensdata,setup)
rb = setup.ross_bytes;
header.nbytes = intcat(ensdata((3:4)+rb)',8) + rb;
header.ndat = ensdata(6+rb);
header.offsets = rb + ...
    intcat(ensdata([rb + 7+2*[0:header.ndat-1]',...
                    rb + 8+2*[0:header.ndat-1]']),8);
header.ids = intcat(ensdata([header.offsets+1,header.offsets+2]),8);


%----------------------------------------------------------
% Fill an adcp data structure with fixed leader information
% from the given ensemble data.
%----------------------------------------------------------
function config = get_fixed_leader(ensdata,setup);

config = struct();

getopt=@(n,opts) opts{n+1};
offset = intcat(ensdata([7:8] + setup.ross_bytes)',8);
dat = double(ensdata(offset + 1 + setup.ross_bytes:end));
%
config.firmware_ver             = str2num([num2str(dat(3)) '.' num2str(dat(4))]);
config.config                   = [dec2bin(dat(6),8) '-' dec2bin(dat(5),8)];
config.beam_angle               = getopt(bitand(dat(6),3),{15,20,30,nan});
config.n_beams                  = getopt(bitand(dat(6),16)==16,{4,5});
config.beam_freq                = getopt(bitand(dat(5),7),...
                                         {75,150,300,600,1200,2400});
config.beam_pattern             = getopt(bitand(dat(5),8)==8,{'concave','convex'});
config.orientation              = getopt(bitand(dat(5),128)==128,{'down','up'});
config.simflag                  = getopt(boolean(dat(7)),{'real','simulated'});
config.n_cells                  = double(dat(10));
config.pings_per_ensemble       = intcat(dat(11:12)',8);
config.cell_size                = intcat(dat(13:14)',8)/100;
config.blank                    = intcat(dat(15:16)',8)/100;
config.prof_mode                = dat(17);
config.corr_threshold           = dat(18);
config.n_codereps               = dat(19);
config.min_pgood                = dat(20);
config.evel_threshold           = intcat(dat(21:22)',8);
config.time_between_ping_groups = sum(double(dat(23:25)) .* [60 1 0.01]');
coord_sys                       = dat(26);
config.coord                    = dec2bin(coord_sys,8);
config.coord_sys                = getopt(bitand(bitshift(coord_sys,-3),3),...
                                         {'beam','instrument','ship','earth'});
config.use_pitchroll            = getopt(bitand(coord_sys,4)==4,{'no','yes'});  
config.use_3beam                = getopt(bitand(coord_sys,2)==2,{'no','yes'});
config.bin_mapping              = getopt(bitand(coord_sys,1)==1,{'no','yes'});
config.xducer_misalign          = intcat(dat(27:28)',8)/100;
config.magnetic_var             = intcat(dat(29:30)',8)/100;
config.sensor_source            = dec2bin(dat(31),8);
config.sensors_avail            = dec2bin(dat(32),8);
config.bin1_dist                = intcat(dat(33:34)',8)/100;
config.xmit_pulse               = intcat(dat(35:36)',8)/100;
config.water_ref_cells          = intcat(dat(37:38)',8);
config.false_target_thres       = dat(39);
config.xmit_lag                 = intcat(dat(41:42)',8);
config.cpu_serialnum            = intcat(dat(43:50)',8);
config.sysbandwidth             = intcat(dat(51:52)',8);
config.syspower                 = dat(53);
config.serialnum                = intcat(dat(55:58)',8);
% I'm not sure if the following will always work...
if isnan(config.beam_angle)
    config.beam_angle               = dat(59);
end
config.ranges = config.bin1_dist + ...
    [0:config.n_cells-1]'*config.cell_size;


%----------------------------------------------------------
% Fill the ith ensemble of an adcp data structure using 
% the given ensemble data.
%----------------------------------------------------------
function [adcp aidx] = get_data(adcp,ensdata,aidx,setup)
header = get_header(ensdata,setup);
config = get_fixed_leader(ensdata,setup);
pr = match_profile(config,adcp);


% Create a new profile if necessary
new_prof = false;
if isempty(adcp);
    adcp = adcp_init(ensdata,setup);
    new_prof = true;
elseif length(adcp) < pr
    adcp(pr) = adcp_init(ensdata,setup);
    new_prof = true;
end

if new_prof
    aidx(end+1) = 0;
    fprintf(['\rNew configuration identified: ' ...
             '%d %.2fm depth cells\n'],...
            adcp(pr).config.n_cells,adcp(pr).config.cell_size);
    print_progress(sum(aidx),setup.nens,1);
end

aidx(pr) = aidx(pr) + 1;
i = aidx(pr);

nc = adcp(pr).config.n_cells;
nb = adcp(pr).config.n_beams;

for j = 2:length(header.ids)
    switch header.ids(j)

      case 128 % variable leader
        % Get ROSS timestamp with variable leader
        if isfield(adcp(pr).config,'ROSS') & adcp(pr).config.ROSS
            adcp(pr).ross_mtime = datenum(...
                double(ensdata(2 + [1:6])) + ...
                [2000 0 0 0 0 ensdata(2+7)]'/100);
        end
        %
        offset = header.offsets(j) + 1;
        vl = double(ensdata(offset:end));
        adcp(pr).number(i)       = intcat(vl(3:4)',8);
        time                 = vl(5:10)';
        time(end)            = time(end) + vl(11)/100;
        adcp(pr).mtime(i)        = datenum(time);
        adcp(pr).number(i)       = adcp(pr).number(i) + 65536*vl(12);
        % adcp(pr).BIT(i)          = intcat(vl(13:14)',8);
        % adcp(pr).ssnd(i)         = intcat(vl(15:16)',8);
        adcp(pr).depth(i)        = intcat(vl(17:18)',8)/10;
        adcp(pr).heading(i)      = intcat(vl(19:20)',8)/100;
        adcp(pr).pitch(i)        = uint2int(intcat(vl(21:22)',8),16)/100;
        adcp(pr).roll(i)         = uint2int(intcat(vl(23:24)',8),16)/100;
        adcp(pr).salinity(i)     = intcat(vl(25:26)',8);
        adcp(pr).temperature(i)  = intcat(vl(27:28)',8)/100;
        adcp(pr).heading_std(i)  = vl(32);
        adcp(pr).pitch_std(i)    = vl(33);
        adcp(pr).roll_std(i)     = vl(34);
        adcp(pr).pressure(i)     = intcat(vl(49:52)',8)/1000;
        % adcp(pr).pressure_var(i) = intcat(vl(53:56)',8)/1000;

      case 256 % velocity
        offset = header.offsets(j) + 2;
        vel = intcat(...
            [ensdata(offset+2*[1:4*nc]-1),...
             ensdata(offset+2*[1:4*nc])],8);
        vel = uint2int(vel,16);
        adcp(pr).vel(:,1:4,i) = reshape(vel,4,nc)'/1000;

      case 512 % correlation
        offset = header.offsets(j) + 2;
        corr = double(ensdata(offset+[1:4*nc]));
        adcp(pr).corr(:,1:4,i) = reshape(corr,4,nc)';

      case 768 % echo intensity
        offset = header.offsets(j) + 2;
        intens = double(ensdata(offset+[1:4*nc]));
        adcp(pr).intens(:,1:4,i) = reshape(intens,4,nc)';

      case 1024 % percent good
        offset = header.offsets(j) + 2;
        pgood = double(ensdata(offset+[1:4*nc]));
        adcp(pr).perc_good(:,1:4,i) = reshape(pgood,4,nc)';

      case 1280 % status
        offset = header.offsets(j) + 2;
        status = double(ensdata(offset+[1:4*nc]));
        adcp(pr).status(:,1:4,i) = reshape(status,4,nc)';

      case 1536 % bottom-track
        offset                 = header.offsets(j) + 1;
        bt                     = double(ensdata(offset:offset+44));
        adcp(pr).bt_range(:,i)     = intcat(...
                                        [bt(17:2:23),bt(18:2:24)],8)/100;
        btvel                      = intcat(...
                                        [bt(25:2:31),bt(26:2:32)],8);
        adcp(pr).bt_vel(:,i)       = uint2int(btvel,16)/1000;
        adcp(pr).bt_corr(:,i)      = bt(33:36);
        adcp(pr).bt_ampl(:,i)      = bt(37:40);
        adcp(pr).bt_perc_good(:,i) = bt(41:44);

      case 2560 % velocity (beam 5)
        offset = header.offsets(j) + 2;
        vel5 = intcat(...
            [ensdata(offset+2*[1:nc]-1),...
             ensdata(offset+2*[1:nc])],8);
        vel5 = uint2int(vel5,16);
        adcp(pr).vel(:,5,i) = vel5/1000;

      case 3072 % echo intensity (beam 5)
        offset = header.offsets(j) + 2;
        int5 = double(ensdata(offset+[1:nc]));
        adcp(pr).intens(:,5,i) = int5;

    end
end

%----------------------------------------------------------
% Check if 2 acdp configurations are the same
%----------------------------------------------------------
function is_same = same_config(c1,c2)
is_same = (length(c1.ranges) == length(c2.ranges)) && ...
          c1.blank == c2.blank && ...
          c1.cell_size == c2.cell_size;


%----------------------------------------------------------
% Match an ADCP configuration to an adcp data structure
%----------------------------------------------------------
function pr = match_profile(config,adcp)
% If adcp is empty, return 1
% If config matches the nth adcp profile, return n
% If config doesn't match a profile, return length(adcp) + 1
for pr = 1:length(adcp)
    if same_config(config,adcp(pr).config)
        return
    end
end
pr = ~isempty(pr) + 1;


%----------------------------------------------------------
% Simple graphical representation of progress
%----------------------------------------------------------
function print_progress(cur,final,force)
n = 50;
str = repmat('-',1,n);
if cur==1 || ...
       mod(cur,floor(final/n))==0 || ...
       cur ==final || ...
       force
    ndone = floor((cur/final)*n);
    str(1:ndone-1) = '=';
    str(max(1,ndone)) = '|';
    str(1) = '[';
    str(end) = ']';
    pct = sprintf('%.0f',100*cur/final);
    str = ['\r' str ' ' pct '%%'];
    fprintf(str);
    if cur==final
        fprintf('\n')
    end
end


%----------------------------------------------------------
% Remove ensembles from an ADCP data structure
%----------------------------------------------------------
function adcp = rm_ens(adcp,rm)
flds = setdiff(fields(adcp),{'files','config'});
for i = 1:length(flds)
    if ndims(adcp.(flds{i})) == 2
        adcp.(flds{i}) = adcp.(flds{i})(:,~rm);
    elseif ndims(adcp.(flds{i})) == 3
        adcp.(flds{i}) = adcp.(flds{i})(:,:,~rm);
    end
end

