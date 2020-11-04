function dxdt = thetaneurons_full(t, x, K, A, e, ameankinverted)  
    dxdt = gather((1 - cos(x)) + (1 + cos(x)).*(e + K * ameankinverted * (A * pulse(x))));
end

