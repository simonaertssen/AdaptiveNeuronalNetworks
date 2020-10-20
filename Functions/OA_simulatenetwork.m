function [TOA, ZOA, b] = OA_simulatenetwork(tnow, tend, IC, p, odeoptions)
    if nargin < 3
        odeoptions = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
    end
    
%     if isvector(IC)
%         OAIC = zeros(1,p.l);
%         for i = 1:p.l
%             OAIC(i) = sum(exp(1i*IC(p.degrees_in == p.k(i)))) / p.P(p.k(i));
%         end
%     elseif isscalar(IC)
%         OAIC = IC*ones(1,p.l);
%     else
%         error('IC might be wrong?')
%     end
%     
    OAIC = zeros(1,p.l);
    for i = 1:p.l
        OAIC(i) = sum(exp(1i*IC(p.degrees_in == p.k(i)))) / p.P(p.k(i));
    end

    [TOA, b] = ode45(@(t,x) MFROA(t,x,p), [tnow, tend], gather(OAIC), odeoptions);
    ZOA = b*p.P(p.k)/p.N;
end

