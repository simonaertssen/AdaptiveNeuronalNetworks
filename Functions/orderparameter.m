function z = orderparameter(x)
% This function computes the complex valued Kuramoto order parameter.
    z = mean(exp(1i*x), 1);
end

