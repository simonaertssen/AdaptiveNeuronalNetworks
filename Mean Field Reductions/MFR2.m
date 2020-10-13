function dzdt = MFR2(t, z, p)
% Here we compute the differential equation for the mean field reduction 
% using the Ott-Antonsen reduction (see Luke 2013).    
    one = -1i*(z-1)*(z-1)*0.5;
    two = (z+1)*(z+1);
    zc = conj(z);
    H = (1 + (z*z + zc*zc)/6 - 4*real(z)/3);
    dzdt = one + 0.5*two*(-p.delta + 1i*p.eta0 + 1i*p.K*H);
end

