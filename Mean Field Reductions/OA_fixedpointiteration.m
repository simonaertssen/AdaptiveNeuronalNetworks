function bhat = OA_fixedpointiteration(b0, p)
    ntimes = 0;
    bhat = b0; bhatold = -b0;

    while norm(bhat - bhatold) > 1.0e-22 && ntimes < 400
        ntimes = ntimes + 1;
        bhatold = bhat;
        
        bhatc = conj(bhat);
        H = (1 + (bhat.*bhat + bhatc.*bhatc)/6 - 4.*real(bhat)/3);
        z = sqrt(-1i*p.delta +    p.eta0 +    p.OA*H);
        bhat = (1 + z)./(1 - z);
        
        if norm(bhat) > 1
            bhat = (1 - z)./(1 + z);
        end
        
    end
%     disp(['Algorithm took ', num2str(ntimes), ' steps'])
end