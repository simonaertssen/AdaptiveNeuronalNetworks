function [tout, xout, A, Kout] = DOPRI_simulatenetwork(ta,tb,x0,h,p,K0) 
    initarray = make_GPUhandle();
    disp("Start simulation.")
    
    % ODE solver parameters:
    npts = round((tb - ta)/h + 1);
    h = (tb - ta)/(npts-1);
    dim = size(x0);
    
    % Network parameters and handles:
    A = initarray(adjacencymatrix(p.degrees_i, p.degrees_o));
    func = @(t, x, K) thetaneurons_full(t, x, K, A, p.e, 1/p.meandegree, p.a_n);
    
    tout = initarray(linspace(ta,tb,npts));
    xout = initarray(zeros(dim(1),npts)); xout(:,1) = x0;
    if nargin > 5
        Kout = initarray(zeros(dim(1),npts)); Kout(:,1) = K0;
    else
        K0 = p.K; Kout = initarray(K0*ones(dim(1),npts)); 
    end
    
    K7 = h*func(ta, x0, K0);
    for i = 1:(npts-1)
        K1 = K7;
        K2 = h*func(tout(i), xout(:,i) + K1*0.2, Kout(i));
        K3 = h*func(tout(i), xout(:,i) + K1*0.075 +  K2*0.225, Kout(i));
        K4 = h*func(tout(i), xout(:,i) + 44*K1/45 - 56*K2/15 + 32*K3/9, Kout(i));
        K5 = h*func(tout(i), xout(:,i) + 19372*K1/6561 - 25360*K2/2187 + 64448*K3/6561 - 212*K4/729, Kout(i));
        K6 = h*func(tout(i), xout(:,i) + 9017*K1/3168 - 355*K2/33 + 46732*K3/5247 + 49*K4/176 - 5103*K5/18656, Kout(i));
        tmp = xout(:,i) + 35*K1/384 + 500*K3/1113 + 125*K4/192 - 2187*K5/6784 + 11*K6/84;
        
        xout(:,i+1) = wrapToPi(tmp);
        K7 = h*func(tout(i), xout(:,i+1), Kout(i));
    end
    disp("Simulation done.")
end

