function J = MFROAJ(z, p)
    zc = conj(z);
    H = (1 + (z.*z + zc.*zc)/6 - 4.*real(z)/3);
    I = -p.delta + 1i*p.eta0 + 1i*p.OA*H;
    J = zeros(p.Mk, p.Mk);

    for r = 1:p.Mk
        for c = 1:p.Mk
            J(r,c) = 0.5*(z(r)+1)^2 * (1i*p.OA(r,c)*(z(c)-2)/3);
            if r == c
                J(r,c) = J(r,c) - 1i*(z(r)-1) + (z(r)+1)*I(r);
            end
        end
    end
    
%     J(1:p.Mk+1:end) = J(1:p.Mk+1:end) + (-1i*(z-1) + (z+1).*I)';
end