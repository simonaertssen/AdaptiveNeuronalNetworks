function [TOA, ZOA] = OA_simulatenetwork(tnow, tend, IC, parameters, odeoptions)
    if nargin < 3
        odeoptions = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
    end
    opts = optimoptions('fsolve', 'Display','off', 'Algorithm', 'Levenberg-Marquardt');
    normP = parameters.P(parameters.k)/parameters.N;
    OAIC = fsolve(@(b) b*normP - IC, randn(1, parameters.l), opts);
    [TOA, b_i] = ode45(@(t,x) MFROA(t,x,parameters), [tnow, tend], OAIC, odeoptions);
    ZOA = b_i*normP;
end

