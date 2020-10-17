function [TOA, ZOA, b] = OA_simulatenetwork(tnow, tend, IC, p, odeoptions)
    if nargin < 3
        odeoptions = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
    end
% %     opts = optimoptions('fsolve', 'Display','off', 'Algorithm', 'Levenberg-Marquardt');
% %     normP = parameters.P(parameters.k)/parameters.N;
% %     size(zeros(1, parameters.l))
% %     size(zeros(1, parameters.l)*normP)
% 
% %     OAIC = fsolve(@(b) b*normP - IC, zeros(1, parameters.l) + 0.01, opts);
% %     size(OAIC)
%     OAIC = zeros(1,p.l);
%     for i = 1:p.l
%         ks = p.k(i)*ones(p.l,1);
%         OAIC(i) = IC(p.degrees_in == p.k(i))*p.P(p.k(i));
%     end
%     normP = p.P(p.k)/p.N;

    OAIC = zeros(1,p.l);
    for i = 1:p.l
%         meanthetaperdegree = IC(p.degrees_in == p.k(i));
    %     OAIC(i) = orderparameter(meanthetaperdegree) / p.P(p.k(i)) * numel(meanthetaperdegree);
        OAIC(i) = sum(exp(1i*IC(p.degrees_in == p.k(i)))) / p.P(p.k(i));
    end

    [TOA, b] = ode45(@(t,x) MFROA(t,x,p), [tnow, tend], OAIC, odeoptions);
    ZOA = b*p.P(p.k)/p.N;
end

