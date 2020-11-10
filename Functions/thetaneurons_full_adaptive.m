function dxdt = thetaneurons_full_adaptive(t, x, K, A, e, ameankinverted)  
    dxdt = gather((1 - cos(x)) + (1 + cos(x)).*(e + K * ameankinverted * (A * pulse(x))));
end

