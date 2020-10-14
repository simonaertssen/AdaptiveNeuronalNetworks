function out = SelkovImplicitHopf(a,b)
    nom = (a + b)*(a + b) * (b/((a+b)*(a+b)))*(b/((a+b)*(a+b)));
    triple = 4*(a + b*b)*(a + b*b)*(a + b*b);
    den = 2*sqrt(triple - (a + a*a - b*b + 2*a*b*b + b*b*b*b)^2);
    out = -1/8 - nom/den;
end

