function thetas = find_ICs(shape, targetIC, sourceIC)
    if nargin < 3
        sourceIC = randn(shape);
    end
    opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt');
    thetas = fsolve(@(theta) targetIC - orderparameter(theta), sourceIC, opts);
    orderparameter(thetas)
end

