function vel = vel_beam2earth(vel,config,q,varargin)

% Define beam 5 weight for towards-instrument velocity
% (1=only beam 5, 0=only beams 1-4)
idx = find(strcmp(varargin,'beam5weight'));
if ~isempty(idx)
    K = varargin{idx+1};
else
    K = 1;      % weight for beam 5
end
K1 = (1-K); % weight for beams 1-4


%% Constants
nb = size(vel,2);
nc = size(vel,1);
nt = size(vel,3);
isConvex = true; %strcmp(a.config.beam_pattern,'convex');
isUp = false; %strcmp(A.config.orientation,'up');

%% Instrument coordinate transformation
% This depends on beam angle and beam configuration
cb = cosd(config.beam_angle);
sb = sind(config.beam_angle);
if isfield(config,'beam2inst')
    % Sometimes we can get this directly from parsed ADCP data
    B2I = [config.beam2inst, zeros(4,config.n_beams-4)];
else
    c = 2*[~(isConvex==isUp), isUp | ~isConvex, ~isConvex, isUp | isConvex] - 1;
    a = 1/(2*sb)*c;
    [a1 a2 a3 a4] = deal(a(1),a(2),a(3),a(4));
    b = 1/(4*cb);
    if nb==5
        %      v1       v2    v3    v4    v5
        B2I = [a1       a2    0     0     0 ; % X
               0        0     a3    a4    0 ; % Y
               K1*b  K1*b     K1*b  K1*b  K ; % Z
               b        b    -b     -b    0]; % E
    elseif nb==4
        B2I = [a1 a2  0   0;
               0  0   a3  a4;
               b  b   b   b
               b  b  -b  -b];
    end
    % Flip towards-transducer velocity for up-facing ADCPs
    if isUp
        B2I(3,:) = -B2I(3,:);
    end
end
% Reshape velocity matrix and transform to instrument coordinates
vb = reshape(permute(vel,[3 1 2]),nc*nt,nb);
vi = (B2I * vb')';

% TODO: Bin mapping

%% Instrument-to-earth transformation
% Create orientation quaternions for reshaped velocity matrix
q_rep = repmat(q(:),nc,1);

% Apply quaternion trasformations to instrument-coordinate data to get
% earth-coordinate data. Preserve 4th dimension (error velocity).
ve = vi;
ve(:,1:3) = rotateframe(q_rep,vi(:,1:3));

% Reshape velocity matrix to original size
vel = permute(reshape(ve',4,nt,nc), [3 1 2]);
