function [tout,xout] = DOPRI_singleneuron(func,ta,tb,x0,h,I)
    disp("Start integration.")
    % ODE solver parameters:
    npts = round((tb - ta)/h + 1);
    h = (tb - ta)/(npts-1);
    dim = size(x0);
    tout = linspace(ta,tb,npts);
    xout = zeros(dim(1),npts); xout(:,1) = x0;
    
    K7 = h*func(ta, x0, I(tout(1)));
    for i = 1:(npts-1)
        K1 = K7;
        K2 = h*func(tout(i), xout(:,i) + K1*0.2, I(tout(i)));
        K3 = h*func(tout(i), xout(:,i) + K1*0.075 +  K2*0.225, I(tout(i)));
        K4 = h*func(tout(i), xout(:,i) + 44*K1/45 - 56*K2/15 + 32*K3/9, I(tout(i)));
        K5 = h*func(tout(i), xout(:,i) + 19372*K1/6561 - 25360*K2/2187 + 64448*K3/6561 - 212*K4/729, I(tout(i)));
        K6 = h*func(tout(i), xout(:,i) + 9017*K1/3168 - 355*K2/33 + 46732*K3/5247 + 49*K4/176 - 5103*K5/18656, I(tout(i)));
        tmp = xout(:,i) + 35*K1/384 + 500*K3/1113 + 125*K4/192 - 2187*K5/6784 + 11*K6/84;
        
        % Unsuccessful bounds
        %tmp(tmp > pi) = tmp(tmp > pi) - 2*pi;
        %tmp(tmp < -pi) = tmp(tmp < -pi) + 2*pi;
        
        xout(:,i+1) = wrapToPi(tmp); % tmp;
        K7 = h*func(tout(i), xout(:,i+1), I(tout(i)));
    end
    disp("End integration.")
end