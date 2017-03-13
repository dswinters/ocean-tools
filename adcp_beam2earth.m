%% adcp_beam2earth.m
% Usage: A = adcp_beam2earth(A)
% Description: Convert velocities in ADCP beam coordinates to earth coordinates
% Inputs: A: ADCP data structure from rdradcp.m
% Outputs: A: ADCP data structure with velocities rotated to earth coordinates
% Author: Dylan Winters
% Created: Mar 12 2017

function A = adcp_beam2earth(A)

% Check that velocities are in beam coordinates
if ~strcmp(A.config.coord_sys,'beam')
    error(['A.config.coord_sys must be ''beam'''...
           ' (currently ''%s'')'],A.config.coord_sys)
end

%% Constants
nbeams = A.config.n_beams;
ncells = A.config.n_cells;
isConvex = strcmp(A.config.beam_pattern,'convex');
isUp = strcmp(A.config.orientation,'up');
hasBT = isfield(A,'bt_vel');
cb = cosd(A.config.beam_angle);
sb = sind(A.config.beam_angle);
h0 = A.config.xducer_misalign;

%% Functions for getting rotation constants
% Internal tilt sensors don't directly measure pitch/roll
st1 =@(t) sind(A.pitch(t)); % sin(tilt1)
ct1 =@(t) cosd(A.pitch(t)); % cos(tilt1)
st2 =@(t) sind(A.roll(t));  % sin(tilt2)
ct2 =@(t) cosd(A.roll(t));  % cos(tilt2)
% This is how we convert tilt1&tilt2 into pitch:
getPitch =@(t) asind(st1(t)*ct2(t) / sqrt(1 - (st1(t)*st2(t))^2));
sp =@(t) sind(getPitch(t)); % sin(pitch)
cp =@(t) cosd(getPitch(t)); % cos(pitch)
sr =@(t) sind(A.roll(t));   % sin(roll)
cr =@(t) cosd(A.roll(t));   % cos(roll)

%% Define some rotation functions
% Add 4th dimension for error velocity
%          X        Y        Z E
rotx=@(d) [1        0        0 0 ; % X
           0 +cosd(d) -sind(d) 0 ; % Y
           0 +sind(d) +cosd(d) 0 ; % Z
           0        0        0 1]; % E

%          X        Y        Z E
roty=@(d) [+cosd(d) 0 +sind(d) 0 ; % X
           0        1        0 0 ; % Y
           -sind(d) 0 +cosd(d) 0 ; % Z
           0        0        0 1]; % E

%          X        Y        Z E
rotz=@(d) [+cosd(d) +sind(d) 0 0 ; % X
           -sind(d) +cosd(d) 0 0 ; % Y
                  0        0 1 0 ; % Z
                  0        0 0 1]; % E

%% Instrument coordinate transformation
c = 2*[isConvex; ~isConvex; isUp&isConvex; ~(isUp&isConvex)]-1;
a = 1/(2*sb);
b = 1/(4*cb);
d = 1/(4*cb);
%          v1     v2      v3     v4
b2i = [c(1)*a c(2)*a       0      0 ; % X
       0           0  c(3)*a c(4)*a ; % Y
       b           b       b      b ; % Z
       d           d      -d     -d]; % E

%% Earth coordinate transformation
i2e =@(t) rotz(A.heading(t) + h0) * rotx(getPitch(t)) * roty(A.roll(t));

%% Bin mapping
% Beam depth scale factors (tilted->flat)
% These depend on pitch, roll, and beam configuration
m12 = @(t) [-sr(t)*cp(t)*[1;1]; 
            sp(t)*[1;1]]; 
m3  = @(t) cp(t).*cr(t);
sd  = @(t) [cb ./ (m3(t).*cb + c.*m12(t).*sb)];

% Get bin numbers for each depth cell per beam:
% (note: this might return invalid bins depending on pitch & roll)
binmap =@(t) fix(sd(t)*[1:ncells]+0.5);
% Extract bin_mapped beam velocities from raw beam velocities:
% (this assumes that the binmap has no invalid bins)
vb_bm  =@(t,bm) [A.east_vel(bm(1,:),t)'  ;
                 A.north_vel(bm(2,:),t)' ;
                 A.vert_vel(bm(3,:),t)'  ;
                 A.error_vel(bm(4,:),t)'];

%% 3-Beam solutions
% 3-beam solutions are done by setting error velocity to 
% zero and solving for the missing beam's velocity.
% From the instrument coordinate transformation matrix we have:
%             E = d(v1 + v2 - v3 - v4)
err = [1;1;-1;-1];
% Solve for the weights to use when beam b is bad:
% (note: this gives weight 1 for beam b, which will be NaN)
weights_3beam =@(b) -sign(err(b))*err;
solve_3beam =@(v,b) nansum(v.*weights_3beam(b));
% Identify deptch cells where this can be done:
use_3beam =@(v) find(sum(isnan(v))==1);
% Identify bad beams in 3-beam-ready depth cells:
nbeam_bad =@(v,cells) nbeams+1 - sum(cumsum(isnan(v(:,cells))));

%% Process
% For each timestep, we will use the rotations defined previously
% to transform data at all depths simultaneously.
East = nan*A.east_vel;
North = nan*A.north_vel;
Vert = nan*A.vert_vel;
Error = nan*A.error_vel;
for t = 1:length(A.mtime)

    %% Get bin-mapped velocities
    bm = binmap(t); % bin-map the depth cells (bin indices)
    rmbin = bm<1 | bm>ncells; % invalid bins
    bm(rmbin) = 1; % (use cell 1 as a placeholder for indexing)
    vb = vb_bm(t,bm); % bin-mapped beam velocity
    vb(rmbin) = NaN; % remove invalid cells

    %% Apply 3-beam solutions where only 1 beam is bad:
    cidx = use_3beam(vb);
    bidx = nbeam_bad(vb,cidx);
    for i = 1:length(cidx)
        vb(bidx(i),cidx(i)) = solve_3beam(vb(:,cidx(i)),bidx(i));
    end

    %% Apply coordinate transformations
    ve = i2e(t)*b2i*vb;
    East(:,t) = ve(1,:);
    North(:,t) = ve(2,:);
    Vert(:,t) = ve(3,:);
    Error(:,t) = ve(4,:);

    if hasBT
        A.bt_vel(:,t) = ...
            1/1000*i2e(t)*b2i*A.bt_vel(:,t);
    end
end

A.east_vel = East;
A.north_vel = North;
A.vert_vel = Vert;
A.error_vel = Error;
A.config.coord_sys = 'earth';
A.config.bin_mapping = 'yes';
A.config.use_3beam = 'yes';
