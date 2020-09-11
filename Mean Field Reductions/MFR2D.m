function dzdt = MFR2D(t, z, p)
% Here we compute the differential equation for the mean field reduction 
% using a 2D vector for meshgrid plotting purposes.
    dzdt = zeros(2,1); x = z(1); y = z(2);
    
    a = x*x - y*y + 2*x + 1;
    b = y*(x+1);
    H = 1 + (x*x - y*y)/3 - 4*x/3;
    prod = p.eta0 + p.K*H;
    
    dzdt(1) = y*(x-1) - b*prod - a*p.delta/2;
    dzdt(2) = -(x*x - y*y - 2*x + 1)/2 + a*prod/2 - b*p.delta;
end

