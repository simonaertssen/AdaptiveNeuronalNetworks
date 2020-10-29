close all; clear all; clc;

%% Find the relations between ICs:
% Setup 
addpath('../Functions');
addpath('../Mean Field Reductions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 13;
export = true;

%% Theta model parameters
tnow = 0; tend = 8;
h = 0.01;

pars.N = 1000;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;
seed = 2; rng(seed);

pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
IC = - pi/2 * ones(pars.N, 1);
p = prepareOAparameters(make_randomparameters(pars, 0.33));

%% 1. From a specific point in the plane Z0 to IC and OAIC:
Z0 = 0.156 - 1i*0.411;

% fsolve would take too long. Take the ICs as from the formula 
% Z0 = exp(1i*x) so x = -i*log(Z0)
findIC = @(length, z) (-1i*log(z)) * ones(length,1);
Z0fromIC = orderparameter(findIC(p.N, Z0))

findOAIC = @(counts, P, z) conj(z * counts ./ P);
Z0fromOAIC = findOAIC(p.kcount, p.P(p.k), Z0)'*p.P(p.k)/p.N

%%
rng(2)
a = rand(100,1) + 1i*rand(100,1);
b = ones(100,1);
a'*b
b'*a

%% 2. From IC to z0 and OAIC
IC = randn(pars.N,1);

% z0 is easy:
Z0 = orderparameter(IC)

% We know this one, gather per degree and 
OAIC = zeros(1,p.l);
for i = 1:p.l
    OAIC(i) = sum(exp(1i*IC(p.degrees_i == p.k(i)))) / (p.P(p.k(i))+1.0e-24);
end
OAIC*p.P(p.k)/p.N

