function dxdt = thetaneurons(t, x, e, KdivN, a)
    p = pulse(x);
    I_sync = sum(p) - p;
    dxdt = (1 - cos(x)) + (1 + cos(x)).*(e + a * KdivN * I_sync);
end

