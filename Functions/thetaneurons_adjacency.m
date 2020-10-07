function dxdt = thetaneurons_adjacency(t, x, e, KdivN, a, A)
    I_sync = A*pulse(x);
    dxdt = (1 - cos(x)) + (1 + cos(x)).*(e + a * KdivN * I_sync);
end

