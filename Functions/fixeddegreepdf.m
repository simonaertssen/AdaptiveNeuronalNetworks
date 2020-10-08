function pdf = fixeddegreepdf(x)
    pdf = dirac(x);
    pdf(pdf == Inf) = 1; % set Inf to finite value
end

