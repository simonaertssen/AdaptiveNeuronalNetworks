function [TOA, ZOA, b] = OA_simulatenetwork(tnow, tend, OAIC, p, odeoptions)
    if nargin < 5 || ~isstruct(odeoptions)
        odeoptions = odeset('RelTol', 1.0e-9,'AbsTol', 1.0e-9);
    end
    
    if ~isfield(odeoptions,'backwards')
        odeoptions.backwards = false;
    end
    
    if isfield(odeoptions,'backwards') && odeoptions.backwards == true
        tback = -1;
        [~, b] = ode45(@(t,x) MFROA(t,x,p), [tnow, tback], gather(OAIC), odeoptions);
        tnow = tback;
        OAIC = b(end,:);
    end

    [TOA, b] = ode45(@(t,x) MFROA(t,x,p), [tnow, tend], gather(OAIC), odeoptions);
    normp = p.P2D(p.k(:,1), p.k(:,2))/p.N;
    if odeoptions.backwards == true
         transientstime = find(TOA >= 0, 1, 'first') - 1;
%          Implement a way to find the OAIC better if we're too far away
%          test = norm(b(transientstime:end, :)*normp - OAIC*normp)
%          if test > 0.01
%          end
         b = b(transientstime:end, :);
    end
    
    ZOA = b*normp;
end

