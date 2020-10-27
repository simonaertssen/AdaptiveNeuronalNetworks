function thetas_0 = find_ICs(shape, targetIC, sourceIC)
    if nargin < 3
        sourceIC = randn(shape);
    end
    opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt', 'MaxIter', 10);
    thetas = fsolve(@(theta) targetIC - orderparameter(theta), sourceIC, opts);
    orderparameter(thetas)
%     [mean, sigma] = fsolve(@(x) targetIC - orderparameter(wrapToPi(x(1) + x(2) .* randn(shape))), [0, 1], opts);
%     
%     function point = evaluate(shape, targetIC, mu, sigma)
%         disp('inside')
%         mu
%         sigma
%         point = norm(targetIC - orderparameter(wrapToPi(mu + sigma .* randn(shape))));
%     end
%         
%     result = fsolve(@(x) evaluate(shape, targetIC, x(1), x(2)), [1, 1], opts);
%     thetas_0 = orderparameter(wrapToPi(result(1) + result(2) .* randn(shape)))
%     targetIC
end

