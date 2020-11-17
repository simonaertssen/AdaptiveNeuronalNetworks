function [tout, xout, K, Kmeans, info] = DOPRI_simulatenetwork_adaptive(ta,tb,x0,h,p,window) 
    initarray = make_GPUhandle();
    disp("Start simulation.")
    
    % ODE solver parameters:
    N = p.N;
    npts = round((tb - ta)/h + 1);
    h = (tb - ta)/(npts-1);
    maxh = 0.1;
    if h > maxh
        error("Window size too small, we're going to miss the 50ms resolutions")
    end
    timestepsback = round(maxh*(1/h));
    
    % Network parameters and handles:
    K = initarray(ones(N, N));
    Kmeans = initarray(zeros(npts,1)); Kmeans(1) = sum(K, 'all')/N + 1.0e-15;
    info = initarray(zeros(npts,1));
    func = @(t, x, K, Kmean) thetaneurons_full_adaptive(t, x, K, p.e, p.a_n, Kmean);
    
    tout = initarray(linspace(ta,tb,npts));
    xout = initarray(zeros(N,npts)); xout(:,1) = x0;
    spiketimes = initarray(zeros(N,timestepsback));
    
    K7 = h*func(ta, x0, K, Kmeans(1));
    for i = 1:(npts-1)
        t = tout(i); x = xout(:,i);
        K1 = K7;
        K2 = h*func(t, x + K1*0.2, K, Kmeans(i));
        K3 = h*func(t, x + K1*0.075 +  K2*0.225, K, Kmeans(i));
        K4 = h*func(t, x + 44*K1/45 - 56*K2/15 + 32*K3/9, K, Kmeans(i));
        K5 = h*func(t, x + 19372*K1/6561 - 25360*K2/2187 + 64448*K3/6561 - 212*K4/729, K, Kmeans(i));
        K6 = h*func(t, x + 9017*K1/3168 - 355*K2/33 + 46732*K3/5247 + 49*K4/176 - 5103*K5/18656, K, Kmeans(i));
        tmp = x + 35*K1/384 + 500*K3/1113 + 125*K4/192 - 2187*K5/6784 + 11*K6/84;
        xout(:,i+1) = wrapToPi(tmp);
        
        mintimeindex = min(i, timestepsback); 
        
        tmp = any(xout(:,i) - xout(:,i+1) > 2*pi - 0.1);
        if tmp == 1
            totallyredundant = 1;
        end
        info(i) = tmp;
        
        spiketimes(xout(:,i) - xout(:,i+1) > 2*pi - 0.1, end) = t;
        
        for tf = 1:N
            for tn = 1:N
                K(tn,tf) = K(tn,tf) + sum(window(spiketimes(tf,1:mintimeindex) - spiketimes(tn,1:mintimeindex)), 'all');
            end
        end
%         K = K + dW;
%         dW(:,:) = 0;
%         if sum(dW, 'all') > 0
%             warning('above 0')
%         end
        Kmeans(i+1) = sum(K, 'all')/N + 1.0e-15;
        
        spiketimes = circshift(spiketimes,-1,2); % Shift circularly to the left
        spiketimes(:, end) = 0; % Clear for the following loop
        
%         This worked
%         lastspiketimes(tmp ~= xout(:,i+1)) = t;
%         K = K + window(lastspiketimes - lastspiketimes');
%         Kmeans(i+1) = sum(K, 'all')/N + 1.0e-15;

        K7 = h*func(t, xout(:,i+1), K, Kmeans(i+1));
    end
    disp("Simulation done.")
end

