%% adcp_orientation.m
% Usage: [h0 p0 r0] = adcp_orientation(v1,v2,v3,v4,spd)
% Description: Attempt to calculate ADCP orientation offsets w.r.t.
%              a moving platform (e.g. a ship).
% Inputs: v (4-by-N) - Near-surface beam velocities (e.g. vertically 
%                   average first 5 bins)
%       spd (1-by-N) - Ship speed vector
%                 R0 - initial head,pitch,roll offset estimate
% Outputs: R0 - head, pitch, and roll offsets
% 
% Author: Dylan Winters
% Created: 2016-09-23


function [R0] = adcp_orientation(v,spd,R0)

nn = ~any(isnan(v)) & ~isnan(spd) & spd >= 1;
v = v(:,nn);
spd = spd(nn);

options = optimset('display','off');
f =@(X) sqrt(sum(sum((v - s2b(spd,X)').^2)));
R0 = fminunc(f,R0,options);

plots = [];
nplots = [2,2];
nplot = 1;

figure('position',[440 325 520 473],'paperpositionmode','auto')

vm = s2b(spd,R0)';
for i = 1:4
    plots(nplot) = subplot(nplots(1),nplots(2),nplot); nplot=nplot+1;    
    plot(spd,v(i,:),'.'), hold on
    plot(spd,vm(i,:),'.')

    if mod(i,2)==1
        ylabel('Beam Vel (m/s)')
    end
    if i>2
        xlabel('Ship Speed (m/s)')
    end

    title(sprintf('Beam %d',i))
    grid on
    xlim([0 3])
    ylim([-3 3])

end


function vm = s2b(spd,X)
% model ADCP beam velocities as a function of ship speed and 
% ADCP orientation 
vm = (spd'*[-1 0 0])*beam_unit(20,X);
end

function A = i2s(X)
% Instrument to ship transformation matrix
    A = rotz(X(1))*rotx(X(2))*roty(X(3));
end

function B = beam_unit(beam_angle,X);
% Beam unit directions (in ship coordinates)
% as a function of beam angle and X, where X =
% [head pitch roll] offsets
%
% +y = ship forward
% +x = ship right
    th = [2 0 1 -1]*90;
    B = nan(3,4);
    u = [0 0 1]';
    for i = 1:4
        Rx = rotx(0);
        Ry = roty(-beam_angle);
        Rz = rotz(th(i));
        Ro = rotz(X(1))*rotx(X(3))*roty(X(2));
        B(:,i) = Ro*Rz*Ry*Rx*u;
    end
end

end

