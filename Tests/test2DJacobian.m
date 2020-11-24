clear all; close all; clc;
% In this script we will test a 2D version of the OA jacobian, to gain
% insight on the stability of fixpoints and to actually find them.

%% Setup:
addpath('../Functions');
addpath('../Mean Field Reductions/');

pars.N = 5000;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
pars.eta0 = 0.4; pars.delta = 0.7; pars.K = 2;
% pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;

seed = 1; rng(seed);
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_randomparameters(pars, 0.33));
