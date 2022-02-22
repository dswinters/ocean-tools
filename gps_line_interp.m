function f_interp = gps_line_interp(gps,dn_source,prefix,field,varargin)
% Interpolate a field from an NMEA prefix without date/time information to
% date/time from a different NMEA prefix. For example, interpolate heading from
% $HEHDT lines to timestamps from $GPRMC lines.
%
% Inputs: gps - gps structure (see nav_read.m)
%         dn_source - NMEA prefix for date/time field (e.g. 'GPRMC')
%         prefix - NMEA prefix for field to be interpolated (e.g. 'HEHDT')
%         field - field to be interpolated (e.g. 'heading')

% First, transform intra-file line numbers into cumulative inter-file line numbers
lnum_max = 0;
prefixes = setdiff(fields(gps),'files');
for i = 1:length(gps.files)
    lnum_max = 0;
    for ip = 1:length(prefixes)
        idx = gps.(prefixes{ip}).fnum == i-1;
        lnum_max = nanmax([lnum_max; max(gps.(prefixes{ip}).lnum(idx))]);
    end
    for ip = 1:length(prefixes)
        idx = gps.(prefixes{ip}).fnum == i;
        gps.(prefixes{ip}).lnum(idx) = gps.(prefixes{ip}).lnum(idx) + lnum_max;
    end
end

% Then we get the line numbers of all lines from both prefixes
n1 = length(gps.(dn_source).lnum);
n2 = length(gps.(prefix).lnum);

% Create an indexing matrix for later use with the sparse() function
idx = [[  1+zeros(1,n1), 2+zeros(1,n2)];            % row 1: field indices (1 or 2)
       [gps.(dn_source).lnum', gps.(prefix).lnum']; % row 2: line numbers
       [1:n1, 1:n2]];                               % row 3: sequential indices
if ~isempty(idx)
    % This is the cool part... Using some sparse matrix magic, efficiently find the
    % preceeding timestamp line for each line containing our desired field.

    % We create a sparse matrix where, for each entry,
    % row = field index from above
    % column = line number from above
    % value = sequential index from above
    s = full(sparse(idx(1,:),idx(2,:),idx(3,:)));
    % We might have something like this:
    % [1 0 0 2 0 0 0 3 ...   first datatype has 1st,2nd,3rd entries on lines 1,4,8
    %  0 1 0 0 0 0 2 0 ...]  second datatype has 1st,2nd entries on lines 2,7

    % The columns with all zeros indicate a line in the GPS file with other data.
    % We can ignore these with the following line:
    s = s(:,~all(s==0));
    % This turns the above example into a new matrix
    % [1 0 2 0 3 ...
    %  0 1 0 2 0...]

    % Delete entries representing lines before our first timestamp line
    dn0 = find(s(1,:)>0,1,'first');
    s = s(:,dn0:end);

    % Shift the second row 
    s = [s(1,:); s(2,2:end) 0];
    % [1 0 2 0 3 ...
    %  1 0 2 0...]

    % Remove zero columns or repeated entries
    s = s(:,~any(s==0));
    % [1 2 ...
    %  1 2 ...]

    % Remove columns representing line nums that don't exist
    s = s(:,[s(1,:)<=length(gps.(dn_source).lnum)] & ...
            [s(2,:)<=length(gps.(prefix).lnum)]);

    dn = gps.(dn_source).dn(s(1,:));
    f = gps.(prefix).(field)(s(2,:));
    [~,uidx] = unique(dn);
    if ismember('angular',varargin)
        f_interp = gps_interp_heading(dn(uidx),f(uidx),gps.(dn_source).dn,'linear');
    else
        f_interp = interp1(dn(uidx),f(uidx),gps.(dn_source).dn);
    end
else
    warning(['Unable to compute timestamps for $' prefix ' data'])
end
