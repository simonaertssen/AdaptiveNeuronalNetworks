function dzdt = MFROA(t, z, p)
% Here we compute the differential equation for the mean field reduction 
% using the formulation for different types of networks, Ott-Antonsen 2017
    one = -1i.*(z-1).*(z-1)/2;
    two = (z+1).*(z+1);
    zc = conj(z);
    H = (1 + (z.*z + zc.*zc)/6 - 4.*real(z)/3);
    HOA = p.OA*H;
    dzdt = one + two.*(-p.delta + 1i*p.eta0 + 1i*HOA)/2;
end

