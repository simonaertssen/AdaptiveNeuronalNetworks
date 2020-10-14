function dxdt = thetaneurons(t, x, e, KdivN, a)
    p = pulse(x);
    dxdt = (1 - cos(x)) + (1 + cos(x)).*(e + a * KdivN * sum(p));
end

