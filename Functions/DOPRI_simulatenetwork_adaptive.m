function [tout, xout, K, Kmeans, info] = DOPRI_simulatenetwork_adaptive(ta,tb,x0,h,p,K0,plastopts) 
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
    
    % Standard integration arrays:
    tout = linspace(ta,tb,npts);
    xout = initarray(zeros(N,npts)); xout(:,1) = x0;
    func = @(t, x, K, Kmean) thetaneurons_full_adaptive(t, x, K, p.e, p.a_n, Kmean);
    
    % Network parameters and handles:
    if nargin > 6
        K = K0;
    else 
        K = initarray(10*ones(N,N)); % For the learning windows
    end
    Kmeans = initarray(zeros(npts,1)); Kmeans(1) = sum(K, 'all')/N + 1.0e-15;
    lastspiketimes = initarray(zeros(N,1));
    eps = 1.0e-15;
    
    % Pasticity options:
    synaptic_plasticity = false;
    synaptic_scaling = false;
    intrnsic_plasticity = false;
    w_i = 0; w_o = 0; KMAX = inf;
    if isstruct(plastopts)
        KMAX = plastopts.KMAX;
        synaptic_plasticity = isstruct(plastopts.SP);
        if isfield(plastopts.SP,'window')
            window = plastopts.SP.window;
            Kupdate = plastopts.SP.Kupdate;
            if isfield(plastopts.SP,'w_i'); w_i = plastopts.SP.w_i; end
            if isfield(plastopts.SP,'w_o'); w_o = plastopts.SP.w_o; end
        end
        synaptic_scaling = isfield(plastopts,'SS');
        intrnsic_plasticity = isfield(plastopts,'IP');
    end
    
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
        
        pulse = xout(:,i) - xout(:,i+1) > 2*pi - 0.1;
        if synaptic_plasticity && any(pulse == 1)
            % Add weight updates if applicable (only Kempter)
            K(:,pulse) = K(:,pulse) + w_i;
            K(pulse,:) = K(pulse,:) + w_o;
            
            lastspiketimes(pulse) = t;
            dW = window(lastspiketimes - lastspiketimes');
            % Filter out all zeros that do not contriubute to the learning:
            % that is where lastspiketimes == 0
            nonzerotimes = find(lastspiketimes)';
            combos = combvec(nonzerotimes,nonzerotimes);
            idx = sub2ind(size(dW),combos(1,:),combos(2,:));
            %K(idx) = K(idx) + dW(idx);
            K(idx) = Kupdate(K(idx),dW(idx));
            if synaptic_scaling
                K = K * Kmeans(i) ./ (sum(K,1) + eps);
            end
            % Rescale to the maximum
            K = min(-KMAX, max(KMAX, K));
        end
        
        if intrnsic_plasticity
            toimplement = 0;
        end
        
        Kmeans(i+1) = (sum(K, 'all') + eps)/N;

        K7 = h*func(t, xout(:,i+1), K, Kmeans(i+1));
    end
    disp("Simulation done.")
end

