function dzdt = MFR2(t, z, p)
% Here we compute the differential equation for the mean field reduction 
% using the Ott-Antonsen reduction (see Luke 2013).

    if isfield(p,'bif') == 1 && t >= p.bif.t
        if isfield(p.bif,'eta0') == 1
            p.e = (p.e - p.eta0)/p.delta * p.bif.delta + p.bif.eta0;
            p.eta0 = p.bif.eta0;
        end

        if isfield(p.bif,'delta') == 1
            p.delta = p.bif.delta;
        end

        if isfield(p,'K') == 1
            p.K = p.bif.K;
        end
    end
    
    one = -1i.*(z-1).*(z-1)/2;
    two = (z+1).*(z+1);
    zc = conj(z);
    H = (1 + (z.*z + zc.*zc)/6 - 4.*real(z)/3);
    dzdt = one + two.*(-p.delta + 1i*p.eta0 + 1i*p.K.*H)/2;
end

