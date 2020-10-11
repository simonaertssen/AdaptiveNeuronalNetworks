function z_pdf = orderparameter_oa(b_i, P, ks, N)
% This function computes the complex valued mean field order parameter as treated in OttAntonsen2017.
    z_pdf = b_i * P(ks)/N;
end

