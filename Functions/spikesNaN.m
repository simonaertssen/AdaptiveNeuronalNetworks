function [tout,xout] = spikesNaN(tout,xout)
    %spikes = bitor(diff(xout, 1, 2) < -2*pi + 0.1, diff(xout, 1, 2) > 2*pi - 0.1);
    spikes = diff(xout, 1, 2) < -2*pi + 0.1;
    xout(spikes) = NaN;
end