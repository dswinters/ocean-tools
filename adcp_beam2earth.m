%% adcp_beam2earth.m
% Usage: A = adcp_beam2earth(A)
% Description: Convert velocities in ADCP beam coordinates to earth coordinates
% Inputs: A: ADCP data structure from rdradcp.m
% Outputs: A: ADCP data structure with velocities rotated to earth coordinates
% Author: Dylan Winters
% Created: Mar 12 2017

function A = adcp_beam2earth(A,varargin)

if ~strcmp(A.config.coord_sys,'beam')
    error(['A.config.coord_sys must be ''beam'''...
           ' (currently ''%s'')'],A.config.coord_sys)
end

% Replace NaN heading with zeros so MATLAB doesn't complain
noheading = isnan(A.heading);
A.heading(noheading) = 0;

% Make a function to pad coordinate transformation matrices
% so we can use matrix products for full transformation
pad = @(m) cat(2,cat(1,m,[0 0 0]),[0;0;0;1]);

% Beam to instrument coordinate transformation matrix
CB = cosd(A.config.beam_angle);
SB = sind(A.config.beam_angle);
c = 2*strcmp(A.config.beam_pattern,'convex')-1; % 1: convex; -1: concave
a = 1/(2*SB);
b = 1/(4*CB);
d = a/sqrt(2);
b2i = [c*a -c*a  0     0  ;
       0    0    -c*a  c*a;
       b    b     b    b  ;
       d    d    -d   -d  ];

% Convert instrument pitch/roll to ship pitch
A.tilt1 = A.pitch;
A.tilt2 = A.roll;
A.pitch = atand(tand(A.tilt1).*cos(A.roll));
A.pitch = asind(sind(A.tilt1).*cosd(A.roll) ./ ...
               sqrt(1 - (sind(A.tilt1).*sind(A.roll)).^2));

% Time-dependent ship-to-earth coordinate transformation matrix
s2e = @(t) pad(rotz(-A.heading(t)));

% Time-dependent instrument-to-ship coordinate transformation matrix
i2s = @(t) pad(                       ...
    rotz(-A.config.xducer_misalign) * ...
    rotx(A.pitch(t))                * ...
    roty(A.roll(t)));

% Time-dependent depth scale factors (for bin mapping)
k = strcmp(A.config.orientation,'up');
c = strcmp(A.config.beam_pattern,'convex');
zs = 2*[c; ~c; k&c; ~(k&c)]-1;
m1 = @(t) [-sind(A.roll(t))*cosd(A.pitch(t))*ones(2,1);
           sind(A.pitch(t))*ones(2,1)];
m2 = @(t)  cosd(A.pitch(t))*cosd(A.roll((t)))*ones(4,1);
ds = @(t) CB ./ (m2(t).*CB + zs.*m1(t).*SB);

% Rotate to earth coordinates
EAST  = nan*A.east_vel;
NORTH = nan*A.east_vel;
VERT  = nan*A.east_vel;
ERR   = nan*A.east_vel;
for t = 1:length(A.mtime)
    % Depth cell mapping
    if A.config.n_cells > 1;
        cells = 1:A.config.n_cells;
        cellidx = fix(feval(ds,t)*cells+0.5);
    else % Sometimes I make a fake adcp structure with 1
         % depth cell to convert BT velocities to earth coords
        cellidx = [1;1;1;1];
    end
    rmcell = cellidx<1 | cellidx>A.config.n_cells;
    cellidx(rmcell) = 1; % placeholder for invalid cells

    % Extract beam velocities from the correct depth cells
    VB = [A.east_vel(cellidx(1,:),t)';
          A.north_vel(cellidx(2,:),t)';
          A.vert_vel(cellidx(3,:),t)';
          A.error_vel(cellidx(4,:),t)'];

    VB(rmcell) = NaN; % remove invalid cells

    % Apply coordinate transformations
    VE = feval(s2e,t)*feval(i2s,t)*b2i*VB;

    % Store earth-coordinate velocities
    EAST(:,t) = VE(1,:);
    NORTH(:,t) = VE(2,:);
    VERT(:,t) = VE(3,:);
    ERROR(:,t) = VE(4,:);
end

A.east_vel = EAST;
A.north_vel = NORTH;
A.vert_vel = VERT;
A.error_vel = ERROR;

% Remove east & north velocities where we have no heading
A.east_vel(:,noheading) = NaN;
A.north_vel(:,noheading) = NaN;

A.config.coord_sys = 'earth';
A.config.bin_mapping = 'yes';