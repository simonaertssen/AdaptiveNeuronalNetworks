function zoa = map_Ztozoa(Z0, counts, P)
    zoa = conj(Z0 * counts ./ P);
end
