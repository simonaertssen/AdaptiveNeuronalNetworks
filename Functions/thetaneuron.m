function dxdt = thetaneuron(t, x, I)
    dxdt = (1 - cos(x)) + (1 + cos(x)).*I;
end

