function [TOA, ZOA, b] = OA_simulatenetwork(tnow, tend, OAIC, p, odeoptions, transients)
    if nargin < 5 || ~isstruct(odeoptions)
        odeoptions = odeset('RelTol', 1.0e-9,'AbsTol', 1.0e-9);
    end
    
    if nargin == 6 && transients == true
        
        tback = -0.5;
        [~, b] = ode45(@(t,x) MFROA(t,x,p), [tnow, tback], gather(OAIC), odeoptions);
        tnow = tback;
        OAIC = b(end,:);
    end

    [TOA, b] = ode45(@(t,x) MFROA(t,x,p), [tnow, tend], gather(OAIC), odeoptions);
    
    if nargin == 6 && transients == true
         transientstime = find(TOA >= 0, 1, 'first') - 1;
         b = b(transientstime:end, :);
    end
    
    ZOA = b*p.P(p.k)/p.N;
end

