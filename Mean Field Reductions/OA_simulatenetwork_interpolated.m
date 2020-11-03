function [TOA, ZOA, b] = OA_simulatenetwork_interpolated(tnow, tend, OAIC, p, transients)
    odeoptions = odeset('RelTol', 1.0e-12,'AbsTol', 1.0e-12);
    
    percentage = 0.5;
    idx = sort(randperm(p.Mk, round(p.Mk*percentage)));
    partly_p = p;
    partly_p.k = p.k(idx);
    partly_p.k_o = p.k_o(idx);
    partly_p.kcount = p.kcount(idx);
    partly_p.OA = p.OA(idx, idx);
    OAIC = OAIC(idx);    
    
%     if nargin == 5 && transients == true
%         tback = -0.5;
%         [~, b] = ode45(@(t,x) MFROA(t,x,p), [tnow, tback], gather(OAIC), odeoptions);
%         tnow = tback;
%         OAIC = b(end,:);
%     end

    [TOA, b] = ode45(@(t,x) MFROA(t,x,partly_p), [tnow, tend], gather(OAIC), odeoptions);
    [tpts, nums] = size(b);
    [Xq, Yq] = meshgrid(1:p.Mk, 1:tpts);
%     size(b')
%     size(Xq)
%     size(Yq)
    X = Xq(:, idx); Y = Yq(:, idx);
%     size(X)
%     size(Y)
    bresized = interp2(X,Y,b,Xq,Yq,'makima');
%     bresized(1)
%     size(bresized)
%     if nargin == 5 && transients == true
%          transientstime = find(TOA >= 0, 1, 'first') - 1;
%          b = b(transientstime:end, :);
%     end
    
    ZOA = bresized*p.P(p.k)/p.N;
end

