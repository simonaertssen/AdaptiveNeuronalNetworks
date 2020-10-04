function z_link = orderparameter_link(x, degrees, A)
% This function computes the complex valued link field order parameter as treated in Timme2017.
    [~, n] = size(x);
    z_link = zeros(1, n);
    degreesum = sum(degrees);
    for i = 1:n
        angles = x(:, i);
        z_link(i) = sum(A.*exp(1i*bsxfun(@minus, angles, angles')), 'all')./degreesum;
    end
end

