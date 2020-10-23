function [TOA, ZOA, b] = OA_simulatenetwork(tnow, tend, IC, p, odeoptions)
    if nargin < 5
        odeoptions = odeset('RelTol', 1.0e-12,'AbsTol', 1.0e-12);
    end
    
    if numel(IC) > 1
        OAIC = zeros(1,p.l);
        for i = 1:p.l
            OAIC(i) = sum(exp(1i*IC(p.degrees_i == p.k(i)))) / (p.P(p.k(i))+1.0e-24);
        end
    elseif numel(IC) == 1
        OAIC = IC*ones(1,p.l);
    else
        error('IC might be wrong?')
    end
    
    [TOA, b] = ode45(@(t,x) MFROA(t,x,p), [tnow, tend], gather(OAIC), odeoptions);
    ZOA = b*p.P(p.k)/p.N;
end

