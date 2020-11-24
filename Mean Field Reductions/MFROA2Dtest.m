function dzdt = MFROA2Dtest(t, z, p)
% Here we compute the differential equation for the mean field reduction 
% using the formulation for different types of networks, Ott-Antonsen 2017,
% in a 2D format (useful for the Jacobian etc).
% Problem: ODE45 only accepts vector format, so watch out with the indexing
    x = real(z); y = imag(z);
    
    frac = 0.5*((x+1).*(x+1) - y.*y);
    KHk = p.OA*(1 + x.*(x - 4)/3);
    
    f = (x-1).*y - frac*p.delta - (x+1).*y.*(p.eta0 + KHk);
    g = -0.5*((x-1).*(x-1) - y.*y) - (x+1).*y*p.delta + frac.*(p.eta0 + KHk);
    dzdt = f + 1i*g;
end

