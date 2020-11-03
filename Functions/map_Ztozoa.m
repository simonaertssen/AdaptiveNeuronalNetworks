function zoa = map_Ztozoa(Z0,p)
    zoa = conj(Z0 * p.kcount ./ (p.P(p.k)));
end
