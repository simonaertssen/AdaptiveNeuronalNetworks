function z0s = drawOAvectors(ICs, cond, p, col)
    [xdim, ydim] = size(ICs);
    z0s = zeros(xdim, ydim);
    tic
    for i = 1:numel(ICs)
        [r, c] = ind2sub([xdim, ydim], i);
        if cond(r,c) == 1
            zoa0 = ones(p.Mk,1)*ICs(i);
            sim1 = DOPRIstep(@(t,x) MFROA(t,x,p),0,zoa0,0.01);
            z0s(r, c) = map_zoatoZ(conj((sim1 - zoa0)'),p);
        end
    end
    q = quiver(real(ICs), imag(ICs), real(z0s), imag(z0s), 0.8, 'color', col);
    q.HandleVisibility = 'off';
end

