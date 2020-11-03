function thetas = find_ICs(sourceIC, targetIC, normalisation)
    opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt', 'MaxIter', 5);
    if nargin < 3 || numel(normalisation) == 1
        thetas = fsolve(@(theta) targetIC - orderparameter(theta), sourceIC, opts);
    else
        thetas = fsolve(@(theta) targetIC - exp(1i*theta)*normalisation, sourceIC, opts);
        thetas = exp(1i*thetas);
    end
end

