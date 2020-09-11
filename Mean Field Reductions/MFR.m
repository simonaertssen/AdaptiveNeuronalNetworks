function dzdt = MFR(t, z, p)
% Here we compute the differential equation for the mean field reduction 
% using the Ott-Antonsen reduction (see Example 2 in Erik Martens 2020). 
% In this formulation an i was lost in the equations, which is found here
% in the last line.

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

    zc = conj(z);
    brackets = (1 + (z.*z + zc.*zc)/6 - 4.*real(z)/3);
    left = p.delta - 1i.*p.eta0 - 1i*p.K.*brackets;
    dzdt = -1/2 * (left .* (1 + z).*(1 + z) + (1 - z).*(1 - z).*1i);
end

