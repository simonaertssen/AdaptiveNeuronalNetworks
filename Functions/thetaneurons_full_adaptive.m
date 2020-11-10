function dxdt = thetaneurons_full_adaptive(t, x, K, e, ameankinverted)  
    dxdt = gather((1 - cos(x)) + (1 + cos(x)).*(e + ameankinverted * (K * pulse(x))));
    
end

