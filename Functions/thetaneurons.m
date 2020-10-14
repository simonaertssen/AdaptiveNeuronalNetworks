function dxdt = thetaneurons(t, x, e, KdivN, a)
    dxdt = (1 - cos(x)) + (1 + cos(x)).*(e + a * KdivN * sum(pulse(x)));
end

