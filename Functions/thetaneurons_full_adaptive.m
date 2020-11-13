function dxdt = thetaneurons_full_adaptive(t, x, K, e, a, Kmean)  
    dxdt = gather((1 - cos(x)) + (1 + cos(x)).*(e + a * (K * pulse(x)))/Kmean);
end

