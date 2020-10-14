clear all; close all; clc;
% Investigate the selkov model on how the Hopf bifurcation arises

%% Setup
pars = struct('a', 0.1, 'b', 0.65);
tnow = 0; tend = 10;
odeoptions = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
[t, x] = ode45(@(t,x) Selkov(t,x,pars), [tnow, tend], 1 + 0.25*1i, odeoptions);

f = figure('Renderer', 'painters', 'Position', [500 1000 400 400]); box on; hold on;
plot(real(x), imag(x))