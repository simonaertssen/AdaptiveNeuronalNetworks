function dW = ChrolCannon2012Window(dt)
    Ap = 0.23; 
    Am = 0.15;
    t_pos = 200;
    t_neg = 2000;
    dt = dt * 1.0e3;
    dW = Ap*exp(-((dt - 15).^2/t_pos)) - Am*exp(-((dt - 20).^2/t_neg));
end
