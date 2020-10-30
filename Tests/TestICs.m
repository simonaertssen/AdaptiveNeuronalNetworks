close all; clear all; clc;

% Terminology: 
% theta0 is the initial condition for the theta model (N real numbers)
% z0 is the initial condition for the OttAntonsen reduced model (Mk complex numbers)
% Z0 is the initial condition in the complex plane (1 complex numbers)

%% Find the relations between initial conditions:
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
theta0 = - pi/2 * ones(pars.N, 1);
p = prepareOAparameters(make_randomparameters(pars, 0.33));

%% 1. From a specific point in the plane Z0 to theta0 and z0:
Z0 = 0.156 - 1i*0.411
assert(numel(Z0) == 1)

% fsolve would take too long. Take the theta0s as from the formula 
% Z0 = exp(1i*x) so x = -i*log(Z0)
findtheta0 = @(length, z) (-1i*log(z)) * ones(length,1);
theta0fromZ0 = findtheta0(p.N, Z0);
Z0fromtheta0 = orderparameter(theta0fromZ0)
assert(numel(theta0fromZ0) == p.N)

findz0 = @(counts, P, z) conj(z * counts ./ P);
z0fromZ0 = findz0(p.kcount, p.P(p.k), Z0);
Z0fromz0 = z0fromZ0'*p.P(p.k)/p.N
assert(numel(z0fromZ0) == p.Mk)

%% 2. From theta0 to Z0 and z0
rng(seed);

theta0 = randn(p.N,1);
assert(numel(theta0) == p.N)

% Z0 is easy:
Z0fromtheta0 = orderparameter(theta0)
assert(numel(Z0fromtheta0) == 1)

% We know how to find z0, gather per degree and multiply per probability
z0fromtheta0 = zeros(1,p.Mk);
for i = 1:p.Mk
    z0fromtheta0(i) = sum(exp(1i*theta0(p.degrees_i == p.k(i)))) / (p.P(p.k(i))+1.0e-24);
end
Z0fromz0 = z0fromtheta0*p.P(p.k)/p.N
assert(numel(z0fromtheta0) == p.Mk)

%% 3. From z0 to Z0 and theta0 
rng(seed);

z0 = randn(p.Mk,1) + randn(p.Mk,1)*1i;

% Z0 is easy:
Z0fromz0 = z0'*p.P(p.k)/p.N

% For theta0 we should repeat the z0 value per degree, randomly distributed
theta0 = zeros(pars.N, 1);
randomindices = randperm(pars.N);
startindex = 1;
testsum = 0;
for i = 1:p.Mk
    indices = p.degrees_i == p.k(i);
    numthetas = p.kcount(i);
    assert(sum(indices) == numthetas)
    testsum = testsum + numthetas; 
    
    endindex = startindex + numthetas;
%     numel(theta0(indices))
%     numel(z0(i)*ones(numthetas,1))
    theta0(indices) = (-1i*log(z0(i)))*ones(numthetas,1);
    startindex = endindex;
end
assert(testsum == p.N)

Z0fromtheta0 = orderparameter(theta0)
