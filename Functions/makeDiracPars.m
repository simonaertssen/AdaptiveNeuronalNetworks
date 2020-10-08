function diracpars = makeDiracPars(pars, netdegree)
    % Make the necessery parameters for the dirac / fixed degree networks
    diracpars = pars;
    diracpars.netdegree = netdegree;
    diracpars.degrees = zeros(pars.N,1);
    degree_idx = randperm(pars.N); 
    diracpars.degrees(randperm(pars.N)) = diracpars.netdegree;

    diracpars.meandegree = diracpars.netdegree;
    diracpars.P = @(x) diracpdf(x - diracpars.netdegree)*pars.N;
end

