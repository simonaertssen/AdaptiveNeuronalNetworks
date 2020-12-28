function phi = Song2017IP(dt)
    T_min = 90; %ms
    T_max = 110; %ms
    lr = 0.012;

    dt = dt * 1.0e3; % Convert to seconds
    phi = zeros(size(dt));
    t_neg_idx = dt < T_min;    
    t_pos_idx = dt > T_max;
    
    phi(t_neg_idx) = -lr*exp((T_min - dt(t_neg_idx))/T_min);
    phi(t_pos_idx) = lr*exp((dt(t_pos_idx) - T_max)/T_max);
end