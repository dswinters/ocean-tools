%% parse_gps.m
%
% Usage
%   gps = parse_gps(files)
%   gps = parse_gps(... , 'progress', uiprogressdlg)
%
% Inputs
% - f_in
%   This can be a filename, cell array of filenames, the output of MATLAB's
%   "dir" command.
%
% Name-value pair arguments
% - 'progress'
%   Specify a UI progress dialog handle to update it with progress information while
%   parsing data.
%
% Outputs
% - gps
%   Data structure containing GPS fields. Each type of NMEA message gets its own
%   sub-structure.
%
% Author
% - Dylan Winters (dylan.winters@oregonstate.edu)

function GPS = parse_gps(f_in,varargin)

    % Make a cell array of filenames from input
    if isstruct(f_in)
        f_in = fullfile({f_in.folder},{f_in.name});
    elseif isstr(f_in)
        f_in = {f_in};
    end

    %% Parse optional inputs
    p = inputParser;
    addParameter(p,'progress',struct(),@(x) isa(x,'matlab.ui.dialog.ProgressDialog'));
    parse(p,varargin{:});
    progress = p.Results.progress;

    % The bulk of the text parsing is done using MATLAB's regexp function.
    % First, create some regexps for types of data we might see in the text.
    % The portions enclosed in parentheses will be extracted, while the rest
    % is just used for matching. The regexp ',?' means that there might be a comma.
    dec  = ',?(-?\d+\.\d+),?'; % positive or negative float w/ decimal places, maybe comma-enclosed
    int  =@(n) sprintf(',?(-?\\d{%d}),?',n); % n-digit positive or negative integer, maybe comma-enclosed
    intu = ',?(\d+),?'; % integer of unknown length, maybe comma-enclosed
    ltr  = ',?[a-zA-Z],?'; % any letter, maybe comma-enclosed
    EW   = ',?([ewEW]),?'; % 'E' or 'W', maybe comma-enclosed

    % Combine the above regexps for single chunks of data into regexps
    % for different types of complete NMEA strings:
    fmt = struct();
    fmt.HEHDT = ['\$HEHDT,' dec];
    fmt.HEROT = ['\$HEROT,' dec];
    fmt.GPGGA = ['\$GPGGA,' ...
                 int(2) int(2) dec ... % time
                 int(2) dec ltr    ... % lat
                 int(3) dec EW     ... % lon
                 intu intu         ... % qual, #sats
                 dec dec           ... % hdop, alt
                 ltr dec];             % alt units, undulation
    fmt.GPHDT = ['\$GPHDT,' dec];
    fmt.GPHEV = ['\$GPHEV,' dec];
    fmt.GPRMC = ['\$GPRMC,' ...
                 int(2) int(2) dec ltr ...
                 int(2) dec ltr        ...
                 int(3) dec EW         ...
                 dec dec               ...
                 int(2) int(2) int(2)];
    fmt.GPZDA = ['\$GPZDA,' ...
                 int(2) int(2) dec ...
                 int(2) int(2) int(4)];
    fmt.PASHR = ['\$PASHR,' ...
                 int(2) int(2) dec ...
                 dec ltr dec dec   ...
                 dec dec dec];
    fmt.PADCP = ['\$PADCP,' intu ...
                 int(4) int(2) int(2) ...
                 int(2) int(2) dec dec];
    fmt.GPVTG = ['\$GPVTG,' dec, ltr, dec, ltr, dec, ltr, dec, ltr];

    %% NMEA prefix-specific substitution filters
    % Replace 'E' and 'W' in GPRMC/GPGGA matrices with '1' and '-1'
    ewflt = struct('str',{'E','W','e','w'},'sub',{'1','-1','1','-1'});
    filts = struct(...
        'GPRMC', ewflt,...
        'GPGGA', ewflt);

    %% Function handles for extracting fields
    %
    % Each file is consecutively parsed for data from each NMEA type. All
    % lines of a single NMEA type are extracted at once, into a matrix D
    % with a row for each line and a column for each raw field. The
    % function handles below provide instructions for converting this
    % matrix into meaningful data.
    %
    % Defining this structure in this way allows for easy looping through NMEA
    % prefixes and fields within each prefix.
    %
    flds = struct(...
        'PASHR',struct(...
            'dn',    @(D) datenum([zeros(size(D,1),3) D(:,1:3)]) ,...
            'head',  @(D) D(:,4)                                 ,...
            'pitch', @(D) D(:,5)                                 ,...
            'roll',  @(D) D(:,6)                                 ,...
            'yaw',   @(D) D(:,7))                                ,...
        'GPGGA',struct(...
            'dn',    @(D) datenum([zeros(size(D,1),3) D(:,1:3)]) ,...
            'lat',   @(D) D(:,4) + D(:,5)/60                     ,...
            'lon',   @(D) D(:,8).*(D(:,6) + D(:,7)/60)           ,...
            'alt',   @(D) D(:,12)                                ,...
            'geoid', @(D) D(:,13))                               ,...
        'HEHDT',struct(...
            'head',  @(D) D(:,1))                                ,...
        'GPHDT',struct(...
            'head',  @(D) D(:,1))                                ,...
        'GPHEV',struct(...
            'heave', @(D) D(:,1))                                ,...
        'HEROT',struct(...
            'rot',   @(D) D(:,1))                                ,...
        'GPRMC',struct(...
            'dn',    @(D) datenum(D(:,[13 12 11 1 2 3])) +        ...
            datenum([2000 0 0 0 0 0])              ,...
            'lat',   @(D) D(:,4) + D(:,5)/60                     ,...
            'lon',   @(D) D(:,8).*(D(:,6) + D(:,7)/60)           ,...
            'speed', @(D) D(:,9) * 0.514444                      ,...
            'course',@(D) D(:,10))                               ,...
        'GPZDA',struct(...
            'dn',    @(D) datenum(D(:,[6 5 4 1 2 3])))           ,...
        'PADCP',struct(...
            'num',   @(D) D(:,1)                                 ,...
            'dn',    @(D) datenum(D(:,[2:7])) - D(:,8)/86400)    ,...
        'GPVTG',struct(...
            'course', @(D) D(:,1)                                 ,...
            'speed',  @(D) D(:,3) * 0.514444));

    nmea_types = fields(flds);

    % % Check the opts struct for manually defined message formats
    % if nargin > 1
    %     if isfield(opts,'fmt')
    %         msgs = fields(opts.fmt);
    %         for i = 1:length(msgs)
    %             fmt.(msgs{i}) = opts.fmt.(msgs{i});
    %         end
    %     end
    %     if isfield(opts,'flds')
    %         msgs = fields(opts.flds);
    %         for i = 1:length(msgs)
    %             flds.(msgs{i}) = opts.flds.(msgs{i});
    %         end
    %     end
    % end

    %% Initialize output structure
    GPS = struct();
    for i = 1:length(nmea_types)
        prefix = nmea_types{i};
        GPS.(prefix) = struct();
        GPS.(prefix).lnum = [];
        GPS.(prefix).fnum = [];
        vars = fields(flds.(prefix));
        for v = 1:length(vars)
            GPS.(prefix).(vars{v}) = [];
        end
    end

    %% Parse!
    for fi = 1:length(f_in)
        [~,fname,~] = fileparts(f_in{fi});
        progress.Message = sprintf('GPS: Processing %s [%d of %d]',fname,fi,length(f_in));
        ftxt = fileread(f_in{fi}); % read entire file text
        for i = 1:length(nmea_types)
            prefix = nmea_types{i};
            %
            [lines, start] = regexp(ftxt,fmt.(prefix),'tokens','start');
            lines = cat(1,lines{:});

            if ~isempty(lines)
                % Apply substitution filters
                if isfield(filts,prefix)
                    for iflt = 1:length(filts.(prefix))
                        lines(strcmp(lines,filts.(prefix)(iflt).str)) = ...
                            {filts.(prefix)(iflt).sub};
                    end
                end
                D = reshape(sscanf(sprintf('%s*',lines{:}),'%f*'),size(lines));
                %
                vars = fields(flds.(prefix));
                % Grab line numbers by counting occurences of newline characters before
                % the start of each line:
                lnum = nan(size(D,1),1);
                lnum(1) = 1 + length(regexp(ftxt(1:start(1)),'\n'));
                for l = 2:length(lnum)
                    lnum(l) = lnum(l-1) + ...
                              length(regexp(ftxt(start(l-1):start(l)),'\n'));
                end
                GPS.(prefix).lnum = cat(1,GPS.(prefix).lnum,lnum);
                GPS.(prefix).fnum = cat(1,GPS.(prefix).fnum,fi*ones(size(D,1),1));
                % Populate struct with variables
                for v = 1:length(vars)
                    GPS.(prefix).(vars{v}) = cat(1,GPS.(prefix).(vars{v}),...
                                                 flds.(prefix).(vars{v})(D));
                end
            end
        end
        progress.Value = fi/length(f_in);
    end

    GPS.files = cell(length(f_in),1);
    for i = 1:length(f_in)
        [~,fname,fext] = fileparts(f_in{i});
        GPS.files{i} = [fname fext];
    end

    % Remove empty fields
    for i = 1:length(nmea_types)
        fld_flds = fields(GPS.(nmea_types{i}));
        has_data = false;
        for j = 1:length(fld_flds)
            has_data = has_data | ~isempty(GPS.(nmea_types{i}).(fld_flds{j}));
        end
        if ~has_data
            GPS = rmfield(GPS,nmea_types{i});
        end
    end

    % Return empty array if no data
    if isempty(setdiff(fields(GPS),{'files'}))
        GPS = [];
    end
end
