function dW = Waddington2014Window(dt)
    lr = 0.1; alpha = 4.0;
    dt = dt * 1.0e3;
    dW = lr*(1 - ((dt-alpha).^2)./alpha^2).*exp(-abs(dt-alpha)./alpha);
end