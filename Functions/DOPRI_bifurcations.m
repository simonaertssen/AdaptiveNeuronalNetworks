function [tout,xout, tdraw,xdraw] = DOPRI_bifurcations(originalfunc,ta,tb,x0,h,p,threshold)
    npts = round((tb - ta)/h + 1);
    h = (tb - ta)/(npts-1);
    dim = size(x0);
    xout = zeros(dim(1),npts); xout(:,1) = x0;
    tout = linspace(ta,tb,npts);
    
    % Make new function handle to improve speed of function evaluation!
    func = @(t, x) originalfunc(t, x, 0, p.e, p.K/p.N, p.a_n);
    
    % Check if parameter bifurcation will be performed:
    bifurcation = false;
    if isfield(p,'bif') == 1
        bifurcation = true;
    end
    
    for i = 1:(npts-1)
        if bifurcation == true && tout(i) >= p.bif.t
            if isfield(p.bif,'eta0') == 1
                p.e = (p.e - p.eta0)/p.delta * p.bif.delta + p.bif.eta0;
                p.eta0 = p.bif.eta0;
            end
            
            if isfield(p.bif,'delta') == 1
                p.delta = p.bif.delta;
            end
            
            if isfield(p,'K') == 1
                p.K = p.bif.K;
            end
            func = @(t, x) originalfunc(t, x, 0, p.e, p.K/p.N, p.a_n);
            bifurcation = false;
        end
        
        K1 = h*func(ta, x0);
        K2 = h*func(tout(i), xout(:,i) + K1/5);
        K3 = h*func(tout(i), xout(:,i) + 3*K1/40 +  9*K2/40);
        K4 = h*func(tout(i), xout(:,i) + 44*K1/45 - 56*K2/15 + 32*K3/9);
        K5 = h*func(tout(i), xout(:,i) + 19372*K1/6561 - 25360*K2/2187 + 64448*K3/6561 - 212*K4/729);
        K6 = h*func(tout(i), xout(:,i) + 9017*K1/3168 - 355*K2/33 + 46732*K3/5247 + 49*K4/176 - 5103*K5/18656);
        tmp = 35*K1/384 + 500*K3/1113 + 125*K4/192 - 2187*K5/6784 + 11*K6/84;
        
        % Register threshold crossing:
        if any(xout(:,i) > threshold)  
            indices = xout(:,i) > pi;
            xout(indices,i) = xout(indices,i) - 2*pi;
        end
        
        xout(:,i+1) = xout(:,i) + tmp;
    end
end