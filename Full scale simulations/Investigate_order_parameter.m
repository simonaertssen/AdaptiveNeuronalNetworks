clear all; close all; clc;
% In this script we will be testing the performance of different order
% parameters as suggested in Timme2017.
% Use a full-scale adjacency matrix simulation for a complete
% investigation.

%% Setup
addpath('../Functions');
addpath('../Mean Field Reductions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 18;

%% Theta model parameters:
tnow = 0; tend = 3;
h = 0.01;

pars.N = 1000;
pars.a_n = 0.666667;
seed = 2; rng(seed);
IC = randn(pars.N, 1);

%% PSR state: one single stable node
pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

% Network distributions and parameters:
fixeddegreepars = make_fixeddegreeparameters(pars, round(pars.N*0.3));
A_fixeddegree = ones(pars.N) - eye(pars.N);
randompars = make_randomparameters(pars, 0.3);
scalefreepars = make_scalefreeparameters(pars, 3);

for i = 1:3
end


% The full scale simulation using the adjacency matrix:
fixeddegreepars = make_fixeddegreeparameters(pars, pars.N - 1);
[t, thetas] = DOPRI_simulatenetwork(tnow,tend,IC,h,fixeddegreepars);
z = orderparameter(thetas);

% The OA mean field theory:
fixeddegreepars = prepareOAparameters(fixeddegreepars);
OAIC = ones(fixeddegreepars.l,1)*z(1) + 0.001*randn(fixeddegreepars.l,1);
[Toa, b_i] = ode45(@(t,x) MFROA(t,x,fixeddegreepars), [tnow, tend], OAIC, options);
Zoa = b_i * fixeddegreepars.P(fixeddegreepars.k)/fixeddegreepars.N;

degrees = sum(A_fixeddegree,2);
z_net = orderparameter_net(thetas, degrees, A_fixeddegree);
z_mf = orderparameter_mf(thetas, degrees);
z_link = orderparameter_link(thetas, degrees, A_fixeddegree);


%%
f_fixeddegree = figure('Renderer', 'painters', 'Position', [50 800 800 200]);

hold on;
plot(t, abs(z), 'LineWidth', 2);
plot(t, abs(z_net), 'LineWidth', 2);
plot(t, abs(z_mf), 'LineWidth', 2);
plot(t, abs(z_link), 'LineWidth', 2);

xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
ylabel('$\| Z (t) \|$','Interpreter','latex', 'FontSize', labelfont)

legend('Kuramoto order parameter', 'Network order parameter', 'Mean field order parameter', 'Link field order parameter', 'FontSize', labelfont-5, 'Location', 'southeast')
removewhitspace();


