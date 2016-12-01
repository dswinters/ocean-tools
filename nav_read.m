%% nav_read.m
% USAGE: 
% NAV = nav_read(f_in,nmea_types)
%
% DESCRIPTION:
% Convert textual NMEA data from 1 or more files into a MATLAB data
% struct.
%
% INPUTS: f_in (cell array), nmea_types (cell array).
% F_IN contains full paths to files containing navigational
% data. NMEA_TYPES contains the prefixes for nmea strings to be read
% (e.g. 'GPGGA', 'GPRMC')
%
% OUTPUTS: NAV (struct).
% NAV contains a field for each type of NMEA string parsed. These
% fields are themselves structs, with subfields for each element
% parsed from the NMEA string. Line numbers and file numbers are also
% saved, with a corresponding list of file names.
% 
% AUTHOR: Dylan Winters [dylan.s.winters@gmail.com]
% CREATED: 2016-09-16

function NAV = nav_read(f_in,nmea_types)

%% Regexps for parsing each type of line
%
% These expressions tell MATLAB how to extract fields from NMEA
% strings.
fmt = struct();
fmt.HEHDT = '\$HEHDT,(\d+\.\d+)';
fmt.HEROT = '\$HEROT,(-?\d+\.\d+)';
fmt.GPGGA = ['\$GPGGA,(\d{2})(\d{2})(\d{2}\.\d+),(\d{2})(\d{2}\.\d+' ...
             '),N,(\d{3})(\d{2}\.\d+),W,\d+,\d+,\d+\.\d+,(-?\d+\.\d+),M,(-?\d+\.\d+).*?'];
fmt.GPHDT = '\$GPHDT,(\d+\.\d+)';
fmt.GPHEV = '\$GPHEV,(-?\d+\.\d+)';
fmt.GPRMC = ['\$GPRMC,(\d{2})(\d{2})(\d{2}\.\d+),A,(\d{2})(\d{2}\.\' ...
             'd+),N,(\d{3})(\d{2}\.\d+),W,(\d+\.\d+),(\d+\.\d+),(\' ...
             'd{2})(\d{2})(\d{2}).*?'];
fmt.GPZDA = '\$GPZDA,(\d{2})(\d{2})(\d{2}\.\d+),(\d{2}),(\d{2}),(\d{4}).*?';
fmt.PASHR = ['\$PASHR,(\d{2})(\d{2})(\d{2}\.\d+),(\d+\.\d+),T,(-?\' ...
             'd+\.\d+),(-?\d+\.\d+),(-?\d+\.\d+),(-?\d+\.\d+),(-?\' ...
             'd+\.\d+).*?'];
fmt.PADCP = '\$PADCP,(\d+),(\d{4})(\d{2})(\d{2}),(\d{2})(\d{2})(\d+\.\d+),(-?\d+\.\d+).*?';

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
        'lon',   @(D) D(:,6) + D(:,7)/60                     ,...
        'alt',   @(D) D(:,8)                                 ,...
        'geoid', @(D) D(:,9))                                ,...
    'HEHDT',struct(...
        'head',  @(D) D(:,1))                                ,...
    'GPHDT',struct(...
        'head',  @(D) D(:,1))                                ,...
    'GPHEV',struct(...
        'heave', @(D) D(:,1))                                ,...
    'HEROT',struct(...
        'rot',   @(D) D(:,1))                                ,...
    'GPRMC',struct(...
        'dn',    @(D) datenum(D(:,[12 11 10 1 2 3]))         ,...
        'lat',   @(D) D(:,4) + D(:,5)/60                     ,...
        'lon',   @(D) D(:,6) + D(:,7)/60)                    ,...
    'GPZDA',struct(...
        'dn',    @(D) datenum(D(:,[6 5 4 1 2 3])))           ,...
    'PADCP',struct(...
        'num',   @(D) D(:,1)                                 ,...
        'dn',    @(D) datenum(D(:,[2:7])) - D(:,8)/86400));

%% Initialize output structure
NAV = struct();
for i = 1:length(nmea_types)
    prefix = nmea_types{i};
    NAV.(prefix) = struct();
    NAV.(prefix).lnum = [];
    NAV.(prefix).fnum = [];
    vars = fields(flds.(prefix));
    for v = 1:length(vars)
        NAV.(prefix).(vars{v}) = [];
    end
end

%% Parse!
for fi = 1:length(f_in)
    disp(sprintf('Processing %s...',f_in{fi}));
    ftxt = fileread(f_in{fi}); % read entire file text
    for i = 1:length(nmea_types)
        prefix = nmea_types{i};
        %
        [lines, start] = regexp(ftxt,fmt.(prefix),'tokens','start');
        lines = cat(1,lines{:});
        if ~isempty(lines)
            D = reshape(sscanf(sprintf('%s*',lines{:}),'%f*'),size(lines));
            %
            vars = fields(flds.(prefix));
            % Grab line numbers first
            lnum = nan(size(D,1),1);
            lnum(1) = 1 + length(regexp(ftxt(1:start(1)),'\n'));
            for l = 2:length(lnum)
                lnum(l) = lnum(l-1) + ...
                          length(regexp(ftxt(start(l-1):start(l)),'\n'));
            end
            NAV.(prefix).lnum = cat(1,NAV.(prefix).lnum,lnum);
            NAV.(prefix).fnum = cat(1,NAV.(prefix).fnum,fi*ones(size(D,1),1));
            % Populate struct with variables
            for v = 1:length(vars)
                NAV.(prefix).(vars{v}) = cat(1,NAV.(prefix).(vars{v}),...
                                             flds.(prefix).(vars{v})(D));
            end
        end
    end
end

