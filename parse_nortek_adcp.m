function adcp = parse_nortek_adcp(dat_or_file,varargin)

%% Parse optional inputs
p = inputParser;
addOptional(p,'parse_cpu_time',false,@(x) islogical(x))
addParameter(p,'progress',struct(),@(x) isa(x,'matlab.ui.dialog.ProgressDialog'));
parse(p,varargin{:});
progress = p.Results.progress;
parse_cpu_time = p.Results.parse_cpu_time;

cpu_time_bytes = 7*parse_cpu_time;
progress = struct();

%% Load binary data
if iscell(dat_or_file)
    % Cell array of filenames
    dat = load_data(dat_or_file,progress);
elseif isstruct(dat_or_file)
    % File info struct from dir()
    dat = load_data(fullfile({dat_or_file.folder},{dat_or_file.name}),progress);
elseif isstring(dat_or_file) | ischar(dat_or_file)
    % Single filename
    dat = load_data({dat_or_file},progress);
else
    % Treat input as raw data
    dat = dat_or_file;
end

%% Find headers and validate their checksums
sync = 0xA5;
family = 0x10;
ids = [0x15 0x16 0x17 0x18 0x1A 0x1B 0x1C 0x1D 0x1E 0x1F 0xA0];
names = strrep({'burst'
                'average'
                'bottom track'
                'interleaved burst'
                'burst altimeter raw'
                'dvl bottom track'
                'echo sounder'
                'dvl water track'
                'altimeter'
                'avg altimeter raw'
                'string'},' ','_');
idmap = containers.Map(num2cell(ids),names);
h = find(dat(1:end-3)==sync & ismember(dat(3:end-1),ids) & dat(4:end)==family);
h = h(h<length(dat)-1);
hlen = double(dat(h+1));
idx = h+hlen<length(dat); % remove incomplete packets
h = h(idx);
hlen = hlen(idx);

id = dat(h+2);
chk = false(length(h),1);
chk_rec = typecast(reshape(dat(h + hlen - 1 + [-1 0])',[],1),'uint16');
chk_comp = uint16(zeros(size(chk_rec)));
for i = 1:length(h)
    [chk(i), chk_comp(i)] = validate_checksum(dat(h(i) + [0:hlen(i)-3]),chk_rec(i));
end
% Discard invalid headers
h = h(chk);
hlen = hlen(chk);
id = id(chk);

uid = unique(id);
% for i = 1:length(uid)
%     disp(sprintf('%s: %d',idmap(uid(i)), sum(id==uid(i))));
% end


%% Validate data checksums
dchk_rec = typecast(reshape(dat(h + hlen - 1 -[3 2])',[],1),'uint16');
dlen = double(typecast(reshape(dat(h + hlen - 1 -[5 4])',[],1),'uint16'));

%% Discard trailing packets
idx = h+hlen+dlen < length(dat);
h = h(idx);
hlen = hlen(idx);
id = id(idx);
dlen = dlen(idx);
dchk_rec = dchk_rec(idx);

dchk = false(size(dchk_rec));
dchksum = 0*dchk;
for i = 1:length(dchk)
    [dchk(i), dchksum(i)] = validate_checksum(dat(h(i) + hlen(i)-1 + [1:dlen(i)]),dchk_rec(i));
end

%% Discard invalid data checksums
h = h(dchk);
hlen = hlen(dchk);
id = id(dchk);
dlen = dlen(dchk);

%% Discard leading/trailing headers
rm = false(size(h));
if parse_cpu_time
    rm = h<cpu_time_bytes+1;
end
rm = rm | (h+hlen+dlen) > length(dat);
h = h(~rm);
hlen = hlen(~rm);
id = id(~rm);
dlen = dlen(~rm);

%% Trim data if number of BT/burst are different
uid = unique(id); % process unique datatypes separately
id_n = ones(size(id));
id_nmax = zeros(size(uid));
for i = 1:length(uid)
    id_n(id==uid(i)) = 1:sum(id==uid(i));
    id_nmax(i) = sum(id==uid(i));
end
n_min = min(id_nmax(ismember(uid,[0x15 0x17])));

kp = id_n <= n_min;
h = h(kp);
hlen = hlen(kp);
id = id(kp);
dlen = dlen(kp);

%% Process data
adcp = struct();
adcp.config = struct();
for i = 1:length(uid)
    % Extract data of this type into a big matrix. Each row is a sample, each
    % column part of some field.
    idx = find(id==uid(i));
    d = @() dat(h(idx) + hlen(idx) + [0:dlen(idx)-1]);

    % Get cpu timestamps if they exist
    if parse_cpu_time
        cpu_time = double(dat(h(idx) + [-cpu_time_bytes:-1]));
        cpu_time(:,6) = cpu_time(:,6) + 1/100*cpu_time(:,7);
        cpu_time(:,1) = cpu_time(:,1) + 2000;
        cpu_time = datenum(cpu_time(:,1:6))';
    end

    % Process data by operating on the columns, depending on the datatype
    switch uid(i)

      %% Burst data record
      case 0x15
        burst = nortek_parse_burst(d());

        % Config fields
        adcp.config.n_beams = burst.nbeams;
        adcp.config.n_cells = burst.ncells;
        adcp.config.depth_cell_length = burst.cell_size(1);
        adcp.config.bin_1_distance = burst.blanking(1);
        adcp.config.serial_number = burst.serial;
        adcp.config.beam_angle = nan;
        adcp.config.frequency = nan;

        % Data fields
        adcp.time = burst.time;
        adcp.cell_depth = adcp.config.bin_1_distance + ...
            [0:adcp.config.n_cells-1]*adcp.config.depth_cell_length;
        adcp.heading = burst.heading;
        adcp.pitch = burst.pitch;
        adcp.roll = burst.roll;
        adcp.vel = burst.vel;
        adcp.echo_intens = burst.amp;
        adcp.corr = burst.cor;
        adcp.ahrs_rotm = burst.ahrs_rotm;
        adcp.ahrs_gyro = burst.ahrs_gyro;
        if parse_cpu_time
            adcp.nuc_time = cpu_time;
        end

      %% Avg data record
      case 0x16
        adcp.nortek_avg = nortek_parse_burst(d());
        if parse_cpu_time
            adcp.nortek_avg.cpu_time = cpu_time;
        end

      %% Bottom track data record
      case 0x17
        bt = nortek_parse_bt(d());
        if parse_cpu_time
            bt.nuc_time = cpu_time;
        end
        adcp.bt_range = bt.dist;
        adcp.bt_vel = bt.vel;
        adcp.bt_time = bt.time;

      case 0x1C
        adcp.ecs = nortek_parse_burst(d());
        if parse_cpu_time
            adcp.ecs.nuc_time = cpu_time;
        end
    end
end
if isfield(adcp,'vel') && isfield(adcp,'bt_vel')
    adcp.vel = adcp.vel - permute(adcp.bt_vel,[3 1 2]);
    adcp.processing.bt_removed_from_vel = true;
end

adcp.config.beam2inst = [1.1831         0   -1.1831         0
                              0   -1.1831         0    1.1831
                         0.5518         0    0.5518         0
                              0    0.5518         0    0.5518];

end

%% Helper functions below
function [chk, chk_comp] = validate_checksum(d,chk_rec)
    len = length(d);
    chk = double(0xB58C); % checksum is initialized to this
                          % sum of all 16-byte segments (remove trailing 8 bytes if odd length)
    chk = chk + sum(double(typecast(d(1:len-mod(len,2)),'uint16')));

    % add the last byte if odd length:
    if mod(len,2)
        chk = chk + bitshift(double(d(end)),8);
    end
    chk_comp = uint16(mod(chk,2^16));
    chk = chk_comp == chk_rec; % final checksum
end

function dat = load_data(files,progress)
%----------------------------------------------------------
% Read and concatenate data in the given list of files
%----------------------------------------------------------
    dat = [];
    for i = 1:length(files)
        if ~exist(files{i},'file')
            error('File not found: %s',files{i})
        end
        [~,fname,fext] = fileparts(files{i});
        progress.Message = sprintf('ADCP: Opening %s %s [%d of %d]',fname,fext,i,length(files));
        fd = fopen(files{i},'r','ieee-le');
        dat = cat(1,dat,uint8(fread(fd,inf,'uint8')));
        fclose(fd);
        progress.Value = i/length(files);
    end
end
