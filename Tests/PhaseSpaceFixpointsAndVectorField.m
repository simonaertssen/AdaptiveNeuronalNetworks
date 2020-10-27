close all; clear all; clc;
addpath('../Functions');
addpath('../Mean Field Reductions/');

%% Test how we can get the fixpoints of the system
clc
pars.N = 10000;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
seed = 2; rng(seed);
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
netp = 0.1;
rdpars = prepareOAparameters(make_randomparameters(pars, netp));

hold on; box on; grid on; axis square;
phasespaceplot();
drawdiraclimitcycle();

IC = find_ICs([pars.N, 1], 0.01 + 1i*0.01);
print(IC)
% eqpts = findeqptsold(wrapToPi(rand(rdpars.l, 1)*2*pi - pi), rdpars.delta, rdpars.eta0, rdpars.K/rdpars.meandegree, rdpars.OA, rdpars.P(rdpars.k)/rdpars.N);
eqpts = findeqptsold(IC, rdpars.delta, rdpars.eta0, rdpars.K/rdpars.meandegree, rdpars.OA, rdpars.P(rdpars.k)/rdpars.N);

% eqpts = findeqpts(wrapToPi(rand(rdpars.l, 1)*2*pi - pi), rdpars)
ZOA = eqpts'*rdpars.P(rdpars.k)/rdpars.N;

scatter(real(ZOA), imag(ZOA), 'r');


function bhat = findeqpts(b0, p)
    ntimes = 0;
    bhat = b0; bhatold = 2*ones(size(b0));
    while norm(bhat - bhatold) > 0.01 && ntimes < 100
        ntimes = ntimes + 1;
        bhatold = bhat;
        orderparameter(bhat)
        
        bhat = MFROA(0, bhat, p);
        
    end
    disp(['Algorithm took ', num2str(ntimes), ' steps'])
end



function bhat = findeqptsold(b0, delta, eta, Kmeank, Passort, v)
    ntimes = 0;
    bhat = b0; bhatold = 2*ones(size(b0));
    while norm(bhat - bhatold) > 0.000001 && ntimes < 100
        ntimes = ntimes + 1;
        bhatold = bhat;
        
        Z = bhat'*v
        plot(real(Z), imag(Z),'r');
        
        z = sqrt(1i * (-delta + 1i*eta + 1i*Kmeank*Passort*bhat));
        bhat = (1 - z)./(1 + z);
        if norm(bhat) < 1; bhat = 1./bhat; end
        
    end
    disp(['Algorithm took ', num2str(ntimes), ' steps'])
end

