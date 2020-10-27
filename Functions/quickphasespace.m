function out = quickphasespace(p)
    addpath('../Functions/');
    addpath('../Mean Field Reductions/');
	p = prepareOAparameters(p);

    figure; hold on; grid on;
    
    % Main grid:
    l = 1; stp = l/2; interval = -l:stp:l;
    [X,Y] = meshgrid(interval,interval);
    th = 0:pi/50:2*pi; 
    in = inpolygon(X, Y, cos(th), sin(th));
    X = X(in); Y = Y(in); %X = X(X > -1 & X < 1); Y = Y(Y > -1 & Y < 1);
    X = reshape(X, numel(X), 1);
    Y = reshape(Y, numel(Y), 1);
    
    odeoptions = odeset('RelTol', 1.0e-2,'AbsTol', 1.0e-5, 'Normcontrol', 'on');
    
    for i = 1:numel(X)
        IC = find_ICs([p.N, 1], X(i) + 1i*Y(i));
        [~, Z] = OA_simulatenetwork(0, 3, IC, p, odeoptions);
        plot(real(Z), imag(Z))    
    end
    scatter(X,Y, 100, 'or', 'filled')
    phasespaceplot();

end

