close all; clear all; clc;
addpath('../Functions');
addpath('../Mean Field Reductions/');

%% Test fixpoint iteration of the logistic map
logistic = @(x) 2.8*x.*(1-x);

xold = 0; xnew = 0.1;
times = 0;
while norm(xnew - xold) > 1.0e-12 && times < 100
    times = times + 1;
    xold = xnew;
    xnew = logistic(xnew);
end

logistic(xnew)
plot(x, logistic(x))
plot(x, x)
scatter(xnew, logistic(xnew))


%% Test how we can get the fixpoints of the system
clc
pars.N = 1000;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;
% pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;

seed = 2; rng(seed);
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
netp = 0.1;
p = prepareOAparameters(make_randomparameters(pars, 0.2));

cla; hold on; box on; grid on; axis square;
phasespaceplot();
% drawdiraclimitcycle();

IC = ones(p.N,1)-randn(pars.N,1);
OAIC = zeros(1,p.l);
for i = 1:p.l
    OAIC(i) = sum(exp(1i*IC(p.degrees_i == p.k(i)))) / (p.P(p.k(i))+1.0e-24);
end

scatter(real(orderparameter(IC)), imag(orderparameter(IC)), 150, '+')
scatter(real(OAIC*p.P(p.k)/p.N), imag(OAIC*p.P(p.k)/p.N), 150, 'x')

eqpts = findeqpts(OAIC', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 'xr');

[~, Z] = OA_simulatenetwork(0, 7, orderparameter(IC), p);
plot(real(Z), imag(Z))

function bhat = findeqpts(b0, p)
    bs = [b0'*p.P(p.k)/p.N];
    ntimes = 0;
    bhat = b0; bhatold = -b0;
    while norm(bhat - bhatold) > 1.0e-22 && ntimes < 200
        ntimes = ntimes + 1;
        bhatold = bhat;
        
        bhatc = conj(bhat);
        H = (1 + (bhat.*bhat + bhatc.*bhatc)/6 - 4.*real(bhat)/3);
%         z = sqrt((-p.delta + 1i*p.eta0 + 1i*p.K*p.OA*H/p.meandegree)/1i);
        z = sqrt((-p.delta + 1i*p.eta0 + 1i*p.K*p.OA*H)/1i);
        bhat = (1 + z)./(1 - z);
        
        disp('norm')
        norm(bhat)
        norm((1 - z)./(1 + z))
        if norm(bhat) > 1
            bhat = (1 - z)./(1 + z);
        end
        
%         bhat = MFROA(0, bhat, p);
%         if norm(bhat) > 1
%             bhat = 1./bhat;
%         end

        bs = [bs, bhat'*p.P(p.k)/p.N];
    end
    disp(['Algorithm took ', num2str(ntimes), ' steps'])
    norm(bhat - bhatold)
    plot(real(bs), imag(bs))
end

