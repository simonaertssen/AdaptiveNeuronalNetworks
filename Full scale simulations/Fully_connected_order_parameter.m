clear all; close all; clc;
% In this script we will be testing the performance of different order
% parameters as suggested in Timme2017.

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

pars.N = 7500;
pars.a_n = 0.666667;
seed = 2; rng(seed);
IC = randn(pars.N, 1);


%% PSR state: one single stable node
pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

% Network distributions and parameters:
diracpars = makeDiracPars(pars, 100);
A_fixeddegree = ones(pars.N) - eye(pars.N);
randompars = makeRandomPars(pars, 0.3);
scalefreepars = makeScalefreePars(pars, 3);

F = @(t, x, e, KdivN, a) thetaneurons_adjacency(t, x, e, KdivN, a, A_fixeddegree);
[t, thetas] = DOPRI_threshold(F, tnow, tend, IC, h, pars);

% Order parameters:
z = orderparameter(thetas);

MRFIC = z(1);
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], MRFIC, options);

figure
hold on;
plot(t, abs(z));
plot(T, abs(Z));


%%
degrees = sum(A_fixeddegree,2);
z_net = orderparameter_net(thetas, degrees, A_fixeddegree);
z_mf = orderparameter_mf(thetas, degrees);
z_link = orderparameter_link(thetas, degrees, A_fixeddegree);


%%
fspace = figure('Renderer', 'painters', 'Position', [50 800 800 200]);

hold on;
plot(t, abs(z), 'LineWidth', 2);
plot(t, abs(z_net), 'LineWidth', 2);
plot(t, abs(z_mf), 'LineWidth', 2);
plot(t, abs(z_link), 'LineWidth', 2);

xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
ylabel('$\| Z (t) \|$','Interpreter','latex', 'FontSize', labelfont)

legend('Kuramoto order parameter', 'Network order parameter', 'Mean field order parameter', 'Link field order parameter', 'FontSize', labelfont-5, 'Location', 'southeast')
removewhitspace();


%% CPW state: the limit cycle
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

% Network distributions and parameters:
diracpars = makeDiracPars(pars, 100);
randompars = makeRandomPars(pars, 0.3);
scalefreepars = makeScalefreePars(pars, 3);

[t, thetas] = DOPRI_threshold(F, tnow, tend, IC, h, pars);

% Order parameters:
z = orderparameter(thetas);
degrees = sum(A_fixeddegree,2);
z_net = orderparameter_net(thetas, degrees, A_fixeddegree);
z_mf = orderparameter_mf(thetas, degrees);
z_link = orderparameter_link(thetas, degrees, A_fixeddegree);

%%
fspace = figure('Renderer', 'painters', 'Position', [50 800 800 200]);

hold on;
plot(t, abs(z), 'LineWidth', 2);
plot(t, abs(z_net), 'LineWidth', 2);
plot(t, abs(z_mf), 'LineWidth', 2);
plot(t, abs(z_link), 'LineWidth', 2);

xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
ylabel('$\| Z (t) \|$','Interpreter','latex', 'FontSize', labelfont)

legend('Kuramoto order parameter', 'Network order parameter', 'Mean field order parameter', 'Link field order parameter', 'FontSize', labelfont-5, 'Location', 'southeast')
removewhitspace();


