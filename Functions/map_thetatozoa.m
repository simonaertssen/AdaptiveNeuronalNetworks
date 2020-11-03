function zoa = map_thetatozoa(theta0,p)
    zoa = zeros(1,p.Mk);
    for i = 1:p.Mk
        zoa(i) = sum(exp(1i*theta0(p.degrees_i == p.k(i)))) / (p.P(p.k(i))+1.0e-24);
    end
end
