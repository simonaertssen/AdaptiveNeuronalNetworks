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

%% 1. From a specific point in the plane Z0 to IC and OAIC:
Z0 = 0.6 - 1i*0.4;

% fsolve would take too long. Take the ICs as from the formula exp(1i*x) = cos(x) + 1i*sin(x)
findIC = @(length, z) ones(length,1) * (acos(real(z)) + 1i*asin(imag(z)));
findIC(10, Z0)
orderparameter(findIC(10, Z0))

