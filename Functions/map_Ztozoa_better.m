function [zoa, IC] = map_Ztozoa_better(Z0,p,method)
    if nargin < 3
        method = 'fmincon';
    end
    Z02D = [real(Z0); imag(Z0)];
    
    % Define the objective function and its bounds:
    objfun = @(z2D) norm(z2D*p.P(p.k)/p.N - Z02D);
    lb = []; ub = lb;
    
    % Define the constraint:
    function [c,ceq] = cnstrnt(x)
        c = x(1,:).^2 + x(2,:).^2 - 1;
        ceq = objfun(x);
    end

    % No other constraints, set those as empty:
    A = []; b = []; Aeq = []; beq = [];

    % Find distribution of Cartesian ICs within the unit circle
    Xs = randn(1, p.Mk)*0.02 + real(Z0);
    Ys = randn(1, p.Mk)*0.02 + imag(Z0);
    magnitude = Xs.^2 + Ys.^2;
    idx = magnitude > 1^2;
    Xs(idx) = Xs(idx) ./ magnitude(idx);
    Ys(idx) = Ys(idx) ./ magnitude(idx);
    % Find the ICs 
    IC = [Xs; Ys];
    
    % Solve:
    if strcmp(method, 'fmincon')
        nonlcon = @cnstrnt;
        fminconopts = optimoptions('fmincon','Display','off','Algorithm','active-set','MaxFunEvals',10*p.Mk,'MaxIter',200);
        zoa = fmincon(@(x)0,IC,A,b,Aeq,beq,lb,ub,nonlcon,fminconopts);
    elseif strcmp(method, 'fsolve')
        fsolveopts = optimoptions('fsolve','Display','iter','Algorithm','levenberg-marquardt','FunValCheck','off');
        zoa = fsolve(objfun,IC,fsolveopts);        
    end
   
    % Return as complex vector
    zoa = zoa(1,:) + 1i*zoa(2,:);
    IC = Xs + 1i*Ys;
    
    check = norm(zoa*p.P(p.k)/p.N - Z0);
    if p.Mk > 800 && check > 0.001
        error("%s found solution with accuracy of %.3e", string(method), check);
    else
        fprintf("%s found solution with accuracy of %.3e\n", string(method), check);
    end
end
