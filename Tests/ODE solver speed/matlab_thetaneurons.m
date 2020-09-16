function dxdt = matlab_thetaneurons(t, x, e, KdivN, a)
    I_sync = sum(matlab_pulse(x));
    dxdt = (1 - cos(x)) + (1 + cos(x)).*(e + a * KdivN * I_sync);
end

