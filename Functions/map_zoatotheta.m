function theta0 = map_zoatotheta(zoa, p)
    theta0 = zeros(p.N, 1);
    for i = 1:p.Mk
        indices = p.degrees_i == p.k(i);
        theta0(indices) = (-1i*log(zoa(i)*p.P(p.k(i))/p.kcount(i)));
    end
end
