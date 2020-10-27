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

hold on
x = 0:0.01:1;
plot(x, logistic(x))
plot(x, x)
scatter(xnew, logistic(xnew))


%% Test how we can get the fixpoints of the system: this works!
clc
pars.N = 1000;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
% pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;
% pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;

seed = 2; rng(seed);
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
netp = 0.1;
p = prepareOAparameters(make_randomparameters(pars, 0.44));

figure; hold on; box on; grid on; axis square;
phasespaceplot();
drawdiraclimitcycle();

% One way:
% IC = -0.8*ones(p.N,1);
% OAIC = zeros(1,p.l);
% for i = 1:p.l
%     OAIC(i) = sum(exp(1i*IC(p.degrees_i == p.k(i)))) / (p.P(p.k(i))+1.0e-24);
% end
% scatter(real(orderparameter(IC)), imag(orderparameter(IC)), 150, '+')
% scatter(real(OAIC*p.P(p.k)/p.N), imag(OAIC*p.P(p.k)/p.N), 150, 'x')
% 
% eqpts = findeqpts(OAIC', p);
% ZOA = eqpts'*p.P(p.k)/p.N;
% scatter(real(ZOA), imag(ZOA), 150, 'xr');

% The other way:
z0 = -0.9 - 1i*0.2;
OAIC = find_ICs(z0*ones(1, p.l), z0, p.P(p.k)/p.N);
scatter(real(z0), imag(z0), 150, '+k')
scatter(real(OAIC*p.P(p.k)/p.N), imag(OAIC*p.P(p.k)/p.N), 150, 'xg');
eqpts = findeqpts(OAIC', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, 'xr');
[~, Z] = OA_simulatenetwork(0, 7, z0, p);
scatter(real(Z(1)), imag(Z(1)), 100, 'ob', 'filled')
plot(real(Z), imag(Z))

z0 = -0.6 - 1i*0.4;
IC = find_ICs(z0*ones(1, p.l), z0, p.P(p.k)/p.N);
scatter(real(z0), imag(z0), 150, '+')
scatter(real(IC*p.P(p.k)/p.N), imag(IC*p.P(p.k)/p.N), 150, 'x');
eqpts = findeqpts(IC', p);
ZOA = eqpts'*p.P(p.k)/p.N;
scatter(real(ZOA), imag(ZOA), 150, 'xr');
[~, Z] = OA_simulatenetwork(0, 7, z0, p);
scatter(real(Z(1)), imag(Z(1)), 100, 'ob', 'filled')
plot(real(Z), imag(Z))

%% Try with fsolve?
z0 = -0.9 - 1i*0.2;
OAIC = find_ICs(z0*ones(1, p.l), z0, p.P(p.k)/p.N);
scatter(real(z0), imag(z0), 150, '+k')
scatter(real(OAIC*p.P(p.k)/p.N), imag(OAIC*p.P(p.k)/p.N), 150, 'xg');

opts = optimset('Display','off', 'Algorithm', 'levenberg-marquardt', 'MaxIter', 10);



%% Functions
function bhat = findeqpts(b0, p)
    bs = b0'*p.P(p.k)/p.N;
    ntimes = 0;
    bhat = b0; bhatold = -b0;

    while norm(bhat - bhatold) > 1.0e-22 && ntimes < 400
        ntimes = ntimes + 1;
        bhatold = bhat;
        
        bhatc = conj(bhat);
        H = (1 + (bhat.*bhat + bhatc.*bhatc)/6 - 4.*real(bhat)/3);
        z = sqrt(-1i*p.delta +    p.eta0 +    p.K*p.OA*H);
        z = sqrt(   (p.delta + 1i*p.eta0 + 1i*p.K*p.OA*H)/1i);
        bhat = (1 + z)./(1 - z);
        
        if norm(bhat) > 1
            bhat = (1 - z)./(1 + z);
        end
        bs = [bs, bhat'*p.P(p.k)/p.N];
        
    end
    disp(['Algorithm took ', num2str(ntimes), ' steps'])
    plot(real(bs), imag(bs))
end