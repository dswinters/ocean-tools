%% parse_adcp.m
%
% Usage
%   adcp = parse_adcp(files)
%   adcp = parse_adcp(files, parse_nuc_timestamps)
%   adcp = parse_adcp(... , 'progress', uiprogressdlg)
%
% Inputs
% - dat_or_file
%   This can be a filename, cell array of filenames, the output of MATLAB's
%   "dir" command, or an array of binary data.
%
% Optional Arguments
% - parse_nuc_timestamps (logical)
% | Flag to parse Nuc timestamps inserted into the ADCP_timestamped* files from
% | ROSE deployments.
%
% Name-value pair arguments
% - 'progress'
%   Specify a UI progress dialog handle to update it with progress information while
%   parsing data.
%
% Outputs
% - adcp
%   Data structure containing ADCP fields.
%
% Author
% - Dylan Winters (dylan.winters@oregonstate.edu)

function adcp = parse_adcp(dat_or_file,varargin)

    %% Parse optional inputs
    p = inputParser;
    addOptional(p,'parse_nuc_timestamps',false,@(x) islogical(x))
    addParameter(p,'progress',struct(),@(x) isa(x,'matlab.ui.dialog.ProgressDialog'));
    parse(p,varargin{:});
    progress = p.Results.progress;
    parse_nuc_timestamps = p.Results.parse_nuc_timestamps;
    nuc_offset = 0;
    if parse_nuc_timestamps
        nuc_offset = 7;
    end

    %% Main parsing routine
    % 1) Load binary data
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

    [h len] = find_headers(dat,progress); % find header locations
    if isempty(h)
        adcp = []; return
    end

    %% Fill adcp data structure with ensemble data
    progress.Message = 'ADCP: Processing data...';
    progress.Indeterminate = 'on';

    % Loop over unique ensemble lengths (there coule be multiple, e.g. for Sentinel
    % V's multi-profile mode)
    [u,~,iu] = unique(len);
    for nl = 1:length(u)
        dat_idx = h(iu==nl) + [-nuc_offset:u(nl)-1];
        % dat_idx = h(iu==nl) + [-nuc_offset:u(nl)+3]; % FIXME: include checksums
        adcp(nl) = parse_ensemble_data(dat(dat_idx));
    end
    progress.Indeterminate = 'off';
    progress.Message = 'ADCP: Processing data... Done!';

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

    function [h len] = find_headers(dat,progress)
    %----------------------------------------------------------
    % Find all header start indices in raw data
    %----------------------------------------------------------
        h = find(dat(1:end-1)==127 & dat(2:end)==127);

        % Need at least 3 bytes after header to get length
        h = h(h < (length(dat)-4));

        % Need at least 7 bytes prior to first header for ROSE timestamp
        if parse_nuc_timestamps
            h = h(h>7);
        end

        if isempty(h)
            h = []; len = []; return
        end

        % Compute ensemble lengths
        len = double(typecast(reshape(dat(h + [2,3])',[],1),'uint16'));
        ensemble_size_max = 5000; % byte limit for ensembles

        % Discard trailing ensembles and ensembles larger than the size limit
        rm = (h + len > length(dat) - 2) | len > ensemble_size_max;
        h = h(~rm);
        len = len(~rm);

        % Verify checksums
        progress.Message = 'ADCP: Verifying ensemble checksums...';
        chk = typecast(reshape(dat(h + len + [0:1])',[],1),'uint16');
        for i = 1:length(chk)
            kp(i) = chk(i) == mod(sum(dat(h(i) + [0:len(i)-1])),65536);
            if mod(i,500)==1
                progress.Value = i/length(chk);
            end
        end
        progress.Value = 1;
        progress.Message = 'ADCP: Verifying ensemble checksums... Done!';

        % Discard header locations with bad checksums
        h = h(kp);
        len = len(kp);

        % Discard packets whose lengths are rare (<10% of data). This should
        % take care any "fake" ensembles in the data caused by coincidental
        % header byte pairs.
        [u,~,iu] = unique(len);
        kp = true(size(h));
        for i = 1:length(u)
            if mean(iu==i) < 0.1
                kp(iu==i) = false;
            end
        end
        h = h(kp);
        len = len(kp);
    end

    function adcp = parse_ensemble_data(dat)
    %----------------------------------------------------------
    % Convert binary data matrix into a structure
    %----------------------------------------------------------
    % Initialize output structure, read number of data types and offsets from first
    % header. If we're parsing bytes inserted prior to each header, account for the
    % number of bytes inserted.
        adcp = struct();
        nuc_offset = 0;
        if parse_nuc_timestamps
            nuc_offset = 7;
        end
        % Read the number of fields (6th byte)
        nfields = dat(1,nuc_offset+6);
        % Compute field data offsets (7th, 9th, ... bytes)
        offsets = nuc_offset + double(typecast(dat(1,nuc_offset+7+[0:2*nfields-1]),'uint16'));
        % The loop below parses fields from each data type.
        for i = 1:nfields

            % For most fields, we can use the following function to extract data
            % colums. It looks nasty, but it's just doing the following:
            % 1) Extract data bytes from all rows (ensembles). The offset plus indices specify colums.
            % 2) Transpose and reshape so we get a long list of bytes, in order.
            % 3) Convert byte sequences into values according to their datatype.
            getdat =@(idx,type) double(typecast(reshape([dat(:,offsets(i)+idx)]',[],1),type))';
            switch dec2hex(typecast(dat(1,offsets(i) + [1:2]),'uint16'),4)
              case '0000' % fixed leader

                % The fixed leader is in every ensemble, but we only need to
                % parse 1. In this case, modify the getdat function to only get
                % data from the first row.
                getdat =@(idx,type) double(typecast(dat(1,offsets(i)+idx),type));

                adcp.config.cpu_fw_ver           = getdat(3,'uint8');
                adcp.config.cpu_fw_rev           = getdat(4,'uint8');
                adcp.config.sys_config           = dec2bin(getdat(5:6,'uint16'),16);
                adcp.config.sym_flag             = logical(getdat(7,'uint8'));
                adcp.config.lag_len              = getdat(8,'uint8');
                adcp.config.n_beams              = getdat(9,'uint8');
                adcp.config.n_cells              = getdat(10,'uint8');
                adcp.config.pings_per_ensemble   = getdat(11:12,'uint16');
                adcp.config.depth_cell_length    = getdat(13:14,'uint16')/100;
                adcp.config.blank_after_transmit = getdat(15:16,'uint16')/100;
                adcp.config.profiling_mode       = getdat(17,'uint8');
                adcp.config.low_corr_thresh      = getdat(18,'uint8');
                adcp.config.n_code_reps          = getdat(19,'uint8');
                adcp.config.perc_good_min        = getdat(20,'uint8');
                adcp.config.error_vel_max        = getdat(21:22,'uint16')/1000;
                adcp.config.tpp_minutes          = getdat(23,'uint8');
                adcp.config.tpp_seconds          = getdat(24,'uint8');
                adcp.config.tpp_hundredths       = getdat(25,'uint8');
                adcp.config.coord_transform      = dec2bin(getdat(26,'uint8'),8);
                adcp.config.heading_alignment    = getdat(27:28,'int16')/100;
                adcp.config.heading_bias         = getdat(29:30,'int16')/100;
                adcp.config.sensor_source        = dec2bin(getdat(31,'uint8'),8);
                adcp.config.sensors_available    = dec2bin(getdat(32,'uint8'),8);
                adcp.config.bin_1_distance       = getdat(33:34,'uint16')/100;
                adcp.config.xmit_puse_length     = getdat(35:36,'uint16')/100;
                adcp.config.wp_ref_layer_avg     = getdat(37:38,'uint8');
                adcp.config.false_target_thresh  = getdat(39,'uint8');
                adcp.config.transmit_lag_dist    = getdat(41:42,'uint16')/100;
                adcp.config.sys_bandwidth        = getdat(51:52,'uint16');
                adcp.config.sys_power            = getdat(53,'uint8');
                adcp.config.serial_number        = getdat(55:58,'uint32');
                adcp.config.beam_angle           = getdat(59,'uint8');

                adcp.config.frequency = 75 * 2^bin2dec(adcp.config.sys_config(end-2:end));
                adcp.cell_depth = adcp.config.bin_1_distance + ...
                    [0:adcp.config.n_cells-1]'*adcp.config.depth_cell_length;


                % The ROSE computers insert timestamps outside of the ADCP's
                % data structure. These don't have an ID or offset, so just grab
                % them along with the fixed leader. Parse inserted bytes after
                % reading fixed leader.
                if parse_nuc_timestamps
                    % This also needs its own getdat in order to specify the data offset manually.
                    getdat =@(idx,type) double(typecast(reshape([dat(:,idx)]',[],1),type))';
                    nuc_time = reshape(getdat(1:7,'uint8'),[nuc_offset, size(dat,1)])';
                    nuc_time(:,1) = nuc_time(:,1) + 2000; % add century
                    nuc_time(:,6) = nuc_time(:,6) + nuc_time(:,7)/100; % add hundredths to seconds
                    adcp.nuc_time = datenum(nuc_time(:,1:6))';
                end

                % FIXME: exploring bad timetsamps
                % nb = adcp.nuc_time < datenum([2020 0 0 0 0 0]);

              case '0080' % variable leader
                ens_num = getdat(3:4,'uint16');
                % clock_year = getdat(5,'uint8');
                % clock_month = getdat(6,'uint8');
                % clock_day = getdat(7,'uint8');
                % clock_hour = getdat(8,'uint8');
                % clock_minute = getdat(9,'uint8');
                % clock_second = getdat(10,'uint8');
                % clock_hundr = getdat(11,'uint8');
                ens_num_msb = getdat(12,'uint8');
                adcp.ens_num = ens_num + 65536 * ens_num_msb;

                % bit_fault = getdat(13,'uint8');
                % bit_reset = getdat(14,'uint8');
                adcp.speed_of_sound = getdat(15:16,'uint16');
                adcp.transducer_depth = getdat(17:18,'uint16')/10;
                adcp.heading = getdat(19:20,'uint16')/100;
                adcp.pitch = getdat(21:22,'int16')/100;
                adcp.roll = getdat(23:24,'int16')/100;
                adcp.salinity = getdat(25:26,'uint16');
                adcp.temperature = getdat(27:28,'int16')/100;
                % mpt_minutes = getdat(29,'uint8');
                % mpt_seconds = getdat(30,'uint8');
                % mpt_hundredths = getdat(31,'uint8');
                adcp.h_std = getdat(32,'uint8')/10;
                adcp.p_std = getdat(33,'uint8')/10;
                adcp.r_std = getdat(34,'uint8')/10;
                % adc_channels = getdat(35:42,'uint8');
                clock_century = getdat(58,'uint8')';
                clock_year = getdat(59,'uint8')';
                clock_month = getdat(60,'uint8')';
                clock_day = getdat(61,'uint8')';
                clock_hour = getdat(62,'uint8')';
                clock_minute = getdat(63,'uint8')';
                clock_second = getdat(64,'uint8')';
                clock_hundr = getdat(65,'uint8')';
                adcp.time = datenum([100*clock_century+clock_year, clock_month, clock_day, ...
                                    clock_hour, clock_minute, clock_second + clock_hundr/100])';

              % These data need some additional reshaping since they have beam,
              % depth, and time dimensions.
              case '0100' % velocity
                adcp.vel = permute(reshape(...
                    getdat(2 + [1:(adcp.config.n_cells * 4 * 2)],'int16')/1000,...
                    4,adcp.config.n_cells,size(dat,1)), [2 1 3]);
              case '0200' % correlation magnitude
                adcp.corr = permute(reshape(...
                    getdat(2 + [1:(adcp.config.n_cells * 4)],'uint8'),...
                    4,adcp.config.n_cells,size(dat,1)), [2 1 3]);
              case '0300' % echo intensity
                adcp.echo_intens = permute(reshape(...
                    getdat(2 + [1:(adcp.config.n_cells * 4)],'uint8'),...
                    4,adcp.config.n_cells,size(dat,1)), [2 1 3]);
              case '0400' % percent good
                adcp.perc_good = permute(reshape(...
                    getdat(2 + [1:(adcp.config.n_cells * 4)],'uint8'),...
                    4,adcp.config.n_cells,size(dat,1)), [2 1 3]);

              % Vertical beam data (Sentinel V)
              case '0A00' % vertical beam velocity
                adcp.vel(:,5,:) = permute(reshape(...
                    getdat(2 + [1:(adcp.config.n_cells * 2)],'int16')/1000,...
                    1,adcp.config.n_cells,size(dat,1)), [2 1 3]);
                adcp.config.n_beams = 5;
              case '0B00' % vertical beam correlation magnitude
                adcp.corr(:,5,:) = permute(reshape(...
                    getdat(2 + [1:(adcp.config.n_cells)],'uint8'),...
                    1,adcp.config.n_cells,size(dat,1)), [2 1 3]);
              case '0C00' % vertical beam echo intensity
                adcp.echo_intens(:,5,:) = permute(reshape(...
                    getdat(2 + [1:(adcp.config.n_cells)],'uint8'),...
                    1,adcp.config.n_cells,size(dat,1)), [2 1 3]);
              case '0D00' % vertical beam percent good
                adcp.perc_good(:,5,:) = permute(reshape(...
                    getdat(2 + [1:(adcp.config.n_cells)],'uint8'),...
                    1,adcp.config.n_cells,size(dat,1)), [2 1 3]);

              % Bottom-track data
              case '0600'
                % BT range is a 3-byte integer with the most significant byte
                % much later in the data packet.
                bt_range_lsb = getdat(17:24,'uint16');
                bt_range_msb = getdat(78:81,'uint8');
                adcp.bt_range = (reshape(bt_range_lsb,4,[]) + 2^16*reshape(bt_range_msb,4,[]))/100;
                adcp.bt_vel = reshape(getdat(25:32,'int16'),4,[])/1000; % mm/s -> m/s
                adcp.bt_perc_good = reshape(getdat(41:44,'uint8'),4,[]);
                adcp.bt_amp = reshape(getdat(41:44,'uint8'),4,[]);
              otherwise
                % disp(dec2hex(typecast(dat(1,offsets(i) + [1:2]),'uint16'),4))
            end
        end
    end % of parse_ensemble_data
end % of parse_adcp
