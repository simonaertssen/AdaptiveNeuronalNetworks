function z_pdf = orderparameter_pdf(x, degrees)
% This function computes the complex valued mean field order parameter as treated in Timme2017.
    z_mf = degrees'*exp(1i*x)./sum(degrees);
end

