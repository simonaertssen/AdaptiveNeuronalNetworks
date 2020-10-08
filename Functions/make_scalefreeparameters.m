function scalefreepars = make_scalefreeparameters(pars, degree, kmin, kmax)
    % Make the necessery parameters for the scalefree networks
    switch nargin
        case 2
            kmin = round(pars.N*3/20);
            kmax = round(pars.N*2/5);
        case 3
            kmax = round(pars.N*2/5);
    end
    if degree < 2
        warning('Scale free networks do not exist for degress less than 2');
    end
    scalefreepars.kmin = kmin;
    scalefreepars.kmax = kmax;
    scalefreepars.degree = degree;
    scalefreepars.P = @(x) scalefreepdf(x, pars.N, scalefreepars.degree, kmin, kmax);

    scalefreepars.degrees_in = randsample(kmin:kmax, pars.N, true, scalefreepars.P(kmin:kmax));

    fsolveoptions = optimset('Display','off');
    scalefreepars.meandegree = fsolve(@(z) scalefreepars.P(z) - mean(scalefreepars.P(kmin:kmax)), kmin, fsolveoptions);
end

