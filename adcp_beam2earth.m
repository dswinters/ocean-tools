%% adcp_beam2earth.m
% Usage: A = adcp_beam2earth(A)
% Description: Convert velocities in ADCP beam coordinates to earth coordinates
% Inputs: A: ADCP data structure from rdradcp.m
% Outputs: A: ADCP data structure with velocities rotated to earth coordinates
% Author: Dylan Winters
% Created: Mar 12 2017

function A = adcp_beam2earth(A)
% clear all, close all
% load('../Data/puget_2017_jan/Spanky/raw/ADCP/SPANKY_2017_01_18_2040.mat')


% Check that velocities are in beam coordinates
if ~strcmp(A.config.coord_sys,'beam')
    error(['A.config.coord_sys must be ''beam'''...
           ' (currently ''%s'')'],A.config.coord_sys)
end

nbeams = A.config.n_beams;
ncells = A.config.n_cells;
isConvex = strcmp(A.config.beam_pattern,'convex');
isUp = strcmp(A.config.orientation,'up');
hasBT = isfield(A,'bt_vel');
cb = cosd(A.config.beam_angle);
sb = sind(A.config.beam_angle);
h0 = A.config.xducer_misalign;
sp =@(t) sind(A.pitch(t));
sr =@(t) sind(A.roll(t));
cp =@(t) cosd(A.pitch(t));
cr =@(t) cosd(A.roll(t));

%% Define some rotation functions
% Add 4th dimension for error velocity
rotx=@(d) [1        0        0 0 ;
           0 +cosd(d) -sind(d) 0 ;
           0 +sind(d) +cosd(d) 0 ;
           0        0        0 1];

roty=@(d) [+cosd(d) 0 +sind(d) 0 ;
           0        1        0 0 ;
           -sind(d) 0 +cosd(d) 0 ;
           0        0        0 1];

rotz=@(d) [+cosd(d) -sind(d) 0 0 ;
           +sind(d) +cosd(d) 0 0 ;
                  0        0 1 0 ;
                  0        0 0 1];

%% Instrument coordinate transformation
c = 2*isConvex-1; % 1: convex; -1: concave
a = 1/(2*sb);
b = 1/(4*cb);
d = a/sqrt(2);
b2i = [c*a -c*a  0     0  ;
       0    0    -c*a  c*a;
       b    b     b    b  ;
       d    d    -d   -d  ];

%% Earth coordinate transformation
% From "adcp coordinate transformation" doc:
% getPitch = atand(tand(A.pitch(t)).*cosd(A.roll(t))); 
% From beam2earth_workhorse.m:
getPitch =@(t) asind(sp(t)*cr(t) / sqrt(1 - (sp(t)*sr(t))^2));
i2s =@(t) rotz(-h0) * rotx(getPitch(t)) * roty(A.roll(t));
s2e =@(t) rotz(-A.heading(t));
% i2e =@(t) rotz(A.heading(t) - h0) * rotx(getPitch(t)) * roty(A.roll(t));

%% Bin mapping
% Beam depth scale factors (tilted->flat)
os = 2*[isConvex; ~isConvex; isUp&isConvex; ~(isUp&isConvex)]-1;
m12 = @(t) [-sr(t)*cp(t)*[1;1]; 
            sp(t)*[1;1]]; 
m3  = @(t) cp(t).*cr(t);
sd  = @(t) [cb ./ (m3(t).*cb + os.*m12(t).*sb)];

% get bin numbers for each depth cell per beam
binmap =@(t) fix(sd(t)*[1:ncells]+0.5);
% extract bin_mapped beam velocities from raw beam velocities
vb_bm  =@(t,bm) [A.east_vel(bm(1,:),t)'  ;
                 A.north_vel(bm(2,:),t)' ;
                 A.vert_vel(bm(3,:),t)'  ;
                 A.error_vel(bm(4,:),t)'];

%% 3-Beam solutions
err = [1;1;-1;-1];
weights_3beam =@(b) -sign(err(b))*err;

East = nan*A.east_vel;
North = nan*A.north_vel;
Vert = nan*A.vert_vel;
Error = nan*A.error_vel;
           
%% Process
for t = 1:length(A.mtime)

    %% Get bin-mapped velocities
    bm = binmap(t); % bin-map the depth cells (bin indices)
    rmbin = bm<1 | bm>ncells; % invalid bins
    bm(rmbin) = 1; % (use cell 1 as a placeholder for indexing)
    vb = vb_bm(t,bm); % bin-mapped beam velocity
    vb(rmbin) = NaN; % remove invalid cells

    % %% Apply 3-beam solutions where only 1 beam is bad:
    % use_3beam = find(sum(isnan(vb))==1);
    % badbeam = nbeams+1 - sum(cumsum(isnan(vb(:,use_3beam))));
    % for i = 1:length(use_3beam)
    %     vb(badbeam(i),use_3beam(i)) = ...
    %         nansum(weights_3beam(badbeam(i)).*vb(:,use_3beam(i)));
    % end

    %% Apply coordinate transformations
    ve = s2e(t)*i2s(t)*b2i*vb;
    % ve = i2e(t)*b2i*vb;
    East(:,t) = ve(1,:);
    North(:,t) = ve(2,:);
    Vert(:,t) = ve(3,:);
    Error(:,t) = ve(4,:);

    if hasBT
        A.bt_vel(:,t) = ...
            1/1000*s2e(t)*i2s(t)*b2i*A.bt_vel(:,t);
            % 1/1000*i2e(t)*b2i*A.bt_vel(:,t);
    end
end

A.east_vel = East;
A.north_vel = North;
A.vert_vel = Vert;
A.error_vel = Error;
A.config.coord_sys = 'earth';
A.config.bin_mapping = 'yes';
A.config.use_3beam = 'yes';
