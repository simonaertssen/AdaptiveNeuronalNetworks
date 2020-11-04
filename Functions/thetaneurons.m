function dxdt = thetaneurons(t, x, e, aKdivN)
    dxdt = gather((1 - cos(x)) + (1 + cos(x)).*(e + aKdivN * sum(pulse(x))));
end

