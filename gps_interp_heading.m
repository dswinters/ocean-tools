function hq = gps_interp_heading(t,h,tq,varargin)
    [t,idx] = unique(t);
    h = h(idx);
    hq = mod(180/pi*angle(interp1(t,cosd(h)+1i*sind(h),tq)),360);
end
