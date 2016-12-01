%% nav_interp_heading.m
% Usage: hq = nav_interp_heading(t,h,tq)
%        hq = nav_interp_heading(t,h,tq,method)
% Description: Interpolate heading in complex space to avoid issues w/ 
%              angle wrapping.
% Inputs: t  - Original time vector
%         h  - Original heading vector
%         tq - Desired time vector
%     method - Use cubic spline interpolation unless another method is specified.
%              currently, only 'linear' is another option.
% Outputs: hq - interpolated heading vector
% 
% Author: Dylan Winters
% Created: 2016-10-21

function hq = nav_interp_heading(t,h,tq,varargin)
    if ismember('linear',varargin)
        hq = mod(180/pi*angle(interp1(t,cosd(h)+1i*sind(h),tq)),360);
    else
        hq = mod(180/pi*angle(spline(t,cosd(h)+1i*sind(h),tq)),360);
    end
end
