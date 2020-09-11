function [tout,xout] = DOPRI_threshold(originalfunc,ta,tb,x0,h,p)
    % ODE solver parameters:
    npts = round((tb - ta)/h + 1);
    h = (tb - ta)/(npts-1);
    dim = size(x0);
    tout = linspace(ta,tb,npts);
    xout = zeros(dim(1),npts); xout(:,1) = x0;
    
    % Make new function handle to improve speed of function evaluation!
    func = @(t, x) originalfunc(t, x, p.e, p.K/p.N, p.a_n);
    
    % Check if we are using the right integration function:
    if isfield(p,'bif') == 1
        disp("Use DOPRI45_bifurcations.m if you wish to study bifurcations")
        return
    end
    
    for i = 1:(npts-1)
        K1 = h*func(ta, x0);
        K2 = h*func(tout(i), xout(:,i) + K1*0.2);
        K3 = h*func(tout(i), xout(:,i) + K1*0.075 +  K2*0.225);
        K4 = h*func(tout(i), xout(:,i) + 44*K1/45 - 56*K2/15 + 32*K3/9);
        K5 = h*func(tout(i), xout(:,i) + 19372*K1/6561 - 25360*K2/2187 + 64448*K3/6561 - 212*K4/729);
        K6 = h*func(tout(i), xout(:,i) + 9017*K1/3168 - 355*K2/33 + 46732*K3/5247 + 49*K4/176 - 5103*K5/18656);
        tmp = xout(:,i) + 35*K1/384 + 500*K3/1113 + 125*K4/192 - 2187*K5/6784 + 11*K6/84;
        
        xout(:,i+1) = wrapToPi(tmp);
    end
end