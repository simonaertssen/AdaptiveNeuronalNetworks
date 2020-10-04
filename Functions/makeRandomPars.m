function randompars = makeRandomPars(pars, netp)
    % Make the necessery parameters for the random networks
    randompars = pars;
    randompars.netp = netp;
    randompars.meandegree = randompars.netp*(pars.N - 1);

    randompars.degrees = poissrnd(randompars.meandegree, [pars.N,1]);
    randompars.P = @(x) poisspdf(x, randompars.meandegree)*pars.N;   
end

