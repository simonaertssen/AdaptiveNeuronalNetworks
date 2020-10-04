function z_net = orderparameter_net(x, A, degreesum)
% This function computes the complex valued net order parameter by Restrepo
% et all, as treated in Timme2017.
    z_i = A*exp(1i*x);
    z_net = sum(z_i)./degreesum;
end

