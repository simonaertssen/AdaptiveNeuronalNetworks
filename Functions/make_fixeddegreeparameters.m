function fixeddegreepars = make_fixeddegreeparameters(pars, netdegree)
    % Make the necessery parameters for the fixeddegree / fixed degree networks
    fixeddegreepars = pars;
    
    if netdegree > pars.N - 1
        error('Net degree too large')
    else 
        fixeddegreepars.netdegree = netdegree;
    end
    fixeddegreepars.degrees_in = zeros(pars.N,1);
    fixeddegreepars.degrees_in(randperm(pars.N)) = fixeddegreepars.netdegree;
    fixeddegreepars.degrees_out = fixeddegreepars.degrees_in(randperm(pars.N));
    
    fixeddegreepars.meandegree = fixeddegreepars.netdegree;
    fixeddegreepars.P = @(x) fixeddegreepdf(x - fixeddegreepars.netdegree)*pars.N;
end

