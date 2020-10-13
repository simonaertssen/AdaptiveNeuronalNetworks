function dxdt = thetaneurons_full(t, x, K, A, e, meankinverted, a)  
    dxdt = (1 - cos(x)) + (1 + cos(x)).*(e + a * K * meankinverted * (A * pulse(x)));
    dxdt = gather(dxdt);
end

