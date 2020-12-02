function dzdt = MFROA2D(t, z, p)
% Here we compute the differential equation for the mean field reduction 
% using the formulation for different types of networks, Ott-Antonsen 2017,
% in a 2D format (useful for the Jacobian etc).
% Problem: ODE45 only accepts vector format, so watch out with the indexing
    x = z(1:2:end); y = z(2:2:end);

    frac = 0.5*((x+1).*(x+1) - y.*y);
    KHk = p.OA*(1 + x.*(x - 4)/3);
    
    dzdt = zeros(2*p.Mk, 1);
    dzdt(1:2:end) = (x-1).*y - frac*p.delta - (x+1).*y.*(p.eta0 + KHk);
    dzdt(2:2:end) = -0.5*((x-1).*(x-1) - y.*y) - (x+1).*y*p.delta + frac.*(p.eta0 + KHk);
end

%     x = z(1,:); y = z(2,:);
% 
%     frac = 0.5*((x+1).*(x+1) - y.*y);
%     KHk = (1 + x.*(x - 4)/3)*p.OA;
%     
%     dzdt = zeros(size(z));
%     dzdt(1,:) = (x-1).*y - frac*p.delta - (x+1).*y.*(p.eta0 + KHk);
%     dzdt(2,:) = -0.5*((x-1).*(x-1) - y.*y) - (x+1).*y*p.delta + frac.*(p.eta0 + KHk);
% end

