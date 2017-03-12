%% adcp_beam2earth.m
% Usage: A = adcp_beam2earth(A)
% Description: Convert velocities in ADCP beam coordinates to earth coordinates
% Inputs: A: ADCP data structure from rdradcp.m
% Outputs: A: ADCP data structure with velocities rotated to earth coordinates
% Author: Dylan Winters
% Created: Mar 12 2017

function A = adcp_beam2earth(A)

% Save beam velocities
if strcmp(A.config.coord_sys,'beam');
    A.beam1_vel = A.east_vel;
    A.beam2_vel = A.north_vel;
    A.beam3_vel = A.vert_vel;
    A.beam4_vel = A.error_vel;
else
    error(['A.config.coord_sys must be ''beam'''...
           ' (currently ''%s'')'],A.config.coord_sys)
end

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
A.pitch_raw = A.pitch;
A.pitch = atand(tand(A.pitch_raw).*cos(A.roll));
A.pitch = asind(sind(A.pitch_raw).*cosd(A.roll) ./ ...
               (1 - (sind(A.pitch_raw).*sind(A.roll)).^2));

% Instrument to ship coordinate transformation matrix
i2s = pad(rotz(A.config.xducer_misalign));

% Time-dependent ship to earth coordinate transformation matrix
s2e = @(t) pad(           ...
    rotz(-A.heading(t)) * ...
    rotx(A.pitch(t))    * ...
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
    % For each ping...

    % Depth cell mapping
    cells = 1:A.config.n_cells;
    cellidx = fix(feval(ds,t)*cells+0.5);
    rmcell = cellidx<1 | cellidx>40;
    cellidx(rmcell) = 1; % placeholder for invalid cells

    % Extract beam velocities from the correct depth cells
    VB = [A.beam1_vel(cellidx(1,:),t)';
          A.beam2_vel(cellidx(2,:),t)';
          A.beam3_vel(cellidx(3,:),t)';
          A.beam4_vel(cellidx(4,:),t)'];
    VB(rmcell) = NaN; % remove invalid cells

    % Apply coordinate transformations
    VE = feval(s2e,t)*i2s*b2i*VB;

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

A.config.coord_sys = 'earth';
A.config.bin_mapping = 'yes';