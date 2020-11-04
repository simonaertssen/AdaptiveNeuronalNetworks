function thetas = find_ICs(sourceIC, targetIC, normalisation)
    function f = fun(x)
        R = targetIC
        f = norm(targetIC - exp(1i*theta)*normalisation);
    end
    
    function [c,ceq] = nonlcon(x)
        c = x(1)^2 + x(2)^2 - 1;
        ceq = [];
    end

    opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt', 'MaxIter', 5);
    if nargin < 3 || numel(normalisation) == 1
        thetas = fsolve(@(theta) targetIC - orderparameter(theta), sourceIC, opts);
    else
%         thetas = fsolve(@(theta) targetIC - exp(1i*theta)*normalisation, sourceIC, opts);
        A = []; b = []; Aeq = []; beq = []; ub = -ones(size(sourceIC)); lb = -ub;
        fun = @(theta) norm(real(targetIC) - exp(1i*theta(1))*normalisation);
        thetas = fmincon(fun,sourceIC,A,b,Aeq,beq,lb,ub);
        thetas = exp(1i*thetas);
    end
end

