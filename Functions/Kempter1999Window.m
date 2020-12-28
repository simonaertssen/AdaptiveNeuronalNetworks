function dW = Kempter1999Window(dt)
    t_syn = 5;
    t_pos = 1;
    t_neg = 20;
    A_p = 1;
    A_n = -1;
    learning_rate = 1.0e-5;
    eps = 1.0e-9;
    dt = -dt * 1.0e3; % Convert to seconds

    dW = zeros(size(dt));
        
    whos dt
    whos dW
    
    t_neg_idx = dt <= 0;
    t_pos_idx = 0 < dt;
    dW(t_neg_idx) = exp(dt(t_neg_idx)/t_syn + eps).*(A_p*(1-dt(t_neg_idx)/t_pos) + A_n*(1-dt(t_neg_idx)/t_neg));
    dW(t_pos_idx) = A_p*exp(-dt(t_pos_idx)/t_pos + eps) + A_n*exp(-dt(t_pos_idx)/t_neg + eps);
    dW = learning_rate * dW;
end