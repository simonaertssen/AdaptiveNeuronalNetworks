function [tout, xout, A, Kout] = DOPRI_simulatenetwork_adaptive(ta,tb,x0,h,p,K0) 
    initarray = make_GPUhandle();
    disp("Start simulation.")
    
    % ODE solver parameters:
    npts = round((tb - ta)/h + 1);
    h = (tb - ta)/(npts-1);
    dim = size(x0);
    
    % Network parameters and handles:
    K = initarray(zeros(p.N, p.N));
    func = @(t, x, K) thetaneurons_full_adaptive(t, x, K, p.e, p.a_n/p.meandegree);
    
    tout = initarray(linspace(ta,tb,npts));
    xout = initarray(zeros(dim(1),npts)); xout(:,1) = x0;
    lastspiketimes = initarray(zeros(dim(1),1));
    lastspiketimeidx = initarray(zeros(dim(1),1, 'logical'));
    
    K7 = h*func(ta, x0, K0);
    for i = 1:(npts-1)
        t = tout(i); x = xout(:,i);
        K1 = K7;
        K2 = h*func(t, x + K1*0.2, K);
        K3 = h*func(t, x + K1*0.075 +  K2*0.225, K);
        K4 = h*func(t, x + 44*K1/45 - 56*K2/15 + 32*K3/9, K);
        K5 = h*func(t, x + 19372*K1/6561 - 25360*K2/2187 + 64448*K3/6561 - 212*K4/729, K);
        K6 = h*func(t, x + 9017*K1/3168 - 355*K2/33 + 46732*K3/5247 + 49*K4/176 - 5103*K5/18656, K);
        tmp = x + 35*K1/384 + 500*K3/1113 + 125*K4/192 - 2187*K5/6784 + 11*K6/84;
        
        xout(:,i+1) = wrapToPi(tmp);
        lastspiketimes(tmp ~= xout(:,i+1)) = t;
        dK = 
        
        K7 = h*func(tout(i), xout(:,i+1), Kout(i));
    end
    disp("Simulation done.")
end

function dW = Waddington2014Window(dt)
    lr = 0.1; alpha = 4.0;
    dW = lr*(1 - ((dt-alpha).^2)./alpha^2).*exp(-abs(dt-alpha)./alpha);
end


function dW = ChrolCannon2012Window(dt)
    Ap = 0.25; Am = 0.1;
    dW = Ap*exp(-((dt - 15).^2/200)) - Am*exp(-((dt - 15).^2/2000));
end
