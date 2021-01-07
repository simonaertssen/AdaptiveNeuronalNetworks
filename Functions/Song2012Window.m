function dW = Song2012Window(dt)
    t_pos = 20;
    t_neg = 20;
    A_pos = 0.1;
    A_neg = -0.12;

    dt = dt * 1.0e3; % Convert to seconds
    dW = zeros(size(dt));
    t_pos_idx = dt > 0;
    t_neg_idx = dt <= 0;
    dW(t_pos_idx) = A_pos*exp(-dt(t_pos_idx)/t_pos);
    dW(t_neg_idx) = A_neg*exp(dt(t_neg_idx)/t_neg);
end