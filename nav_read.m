%% nav_read.m
% USAGE: 
% NAV = nav_read(f_in,nmea_types)
%
% DESCRIPTION:
% Convert textual NMEA data from 1 or more files into a MATLAB data
% struct.
%
% EXAMPLES:
%   1) Read GPGGA data from file1.n1r:
%      nav = nav_read({'file1.n1r'},{'GPGGA'})
%   2) Read GPGGA and HEHDT data from file1.n1r and file2.n1r
%      nav = nav_read({'file1.n1r','file2.n1r'},{'GPGGA','HEHDT'})
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
             int(2) int(2) dec ...
             int(2) dec ltr    ...
             int(3) dec EW     ...
             intu intu dec dec ltr];
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

%% NMEA prefix-specific substitution filters
% Replace 'E' and 'W' in GPRMC/GPGGA matrices with '1' and '-1'
ewflt = struct('str',{'E','W'},'sub',{'1','-1'});
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
        'alt',   @(D) D(:,11)                                ,...
        'geoid', @(D) D(:,12))                               ,...
    'HEHDT',struct(...
        'head',  @(D) D(:,1))                                ,...
    'GPHDT',struct(...
        'head',  @(D) D(:,1))                                ,...
    'GPHEV',struct(...
        'heave', @(D) D(:,1))                                ,...
    'HEROT',struct(...
        'rot',   @(D) D(:,1))                                ,...
    'GPRMC',struct(...
        'dn',    @(D) datenum(D(:,[13 12 11 1 2 3])) + ...
                      datenum([2000 0 0 0 0 0])              ,...
        'lat',   @(D) D(:,4) + D(:,5)/60                     ,...
        'lon',   @(D) D(:,8).*(D(:,6) + D(:,7)/60))          ,...
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

        % Apply substitution filters
        if isfield(filts,prefix)
            for iflt = 1:length(filts.(prefix))
                lines(strcmp(lines,filts.(prefix)(iflt).str)) = ...
                    {filts.(prefix)(iflt).sub};
            end
        end

        if ~isempty(lines)
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

NAV.files = cell(length(f_in),1);
for i = 1:length(f_in)
    [~,fname,fext] = fileparts(f_in{i});
    NAV.files{i} = [fname fext];
end