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
Z0 = 0.01 - 1i*0.4;

% fsolve would take too long. Take the ICs as from the formula 
% Z0 = exp(1i*x) so x = -i*log(Z0)
findIC = @(length, z) (-1i*log(z)) * ones(length,1);
orderparameter(findIC(100, Z0))
% This zorks for IC and OAIC

%% 2. From IC to z0 and OAIC
IC = randn(pars.N,1);

% z0 is easy:
z0 = orderparameter(IC)

% We know this one, gather per degree and 


