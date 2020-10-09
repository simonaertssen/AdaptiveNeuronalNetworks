function dxdt = thetaneurons_full(t, x, K, A, e, Ninverted, a)
    I_sync = A*pulse(x);
    dxdt = (1 - cos(x)) + (1 + cos(x)).*(e + a * K .* I_sync * Ninverted);
end

