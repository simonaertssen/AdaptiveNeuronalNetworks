function fixeddegreepars = make_fixeddegreeparameters(pars, netdegree)
    % Make the necessery parameters for the fixeddegree / fixed degree networks
    fixeddegreepars = pars;
    fixeddegreepars.netdegree = netdegree;
    fixeddegreepars.degrees = zeros(pars.N,1);
    fixeddegreepars.degrees(randperm(pars.N)) = fixeddegreepars.netdegree;

    fixeddegreepars.meandegree = fixeddegreepars.netdegree;
    fixeddegreepars.P = @(x) fixeddegreepdf(x - fixeddegreepars.netdegree)*pars.N;
end

