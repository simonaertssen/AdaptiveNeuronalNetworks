function scalefreepars = makeScalefreePars(pars, degree, kmin, kmax)
    % Make the necessery parameters for the scalefree networks
    switch nargin
        case 2
            kmin = round(pars.N*3/20);
            kmax = round(pars.N*2/5);
        case 3
            kmax = round(pars.N*2/5);
    end

    scalefreepars.degree = degree;
    scalefreepars.P = @(x) scalefreepdf(x, pars.N, scalefreepars.degree, kmin, kmax);

    x = kmin:kmax;
    scalefreepars.degrees = kmin + scalefreepars.P(x)';
    scalefreepars.degrees = randi([kmin, kmax], pars.N, 1);
    
    fsolveoptions = optimset('Display','off');
    scalefreepars.meandegree = fsolve(@(z) scalefreepars.P(z) - mean(scalefreepars.P(kmin:kmax)), kmin, fsolveoptions);
end

