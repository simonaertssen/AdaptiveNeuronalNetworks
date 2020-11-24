function [x, A] = scalefreepdf(x, N, exponent, kmin, kmax)
    xposidx = (x >= kmin & x <= kmax);
    xpos = x(xposidx);
    pdfxpos = xpos.^(-exponent);
    
    normalization = kmin:kmax;
    A = N/sum(normalization.^(-exponent));

    x(xposidx) = A*pdfxpos;
    x(~xposidx) = 0;
end
