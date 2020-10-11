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
labelfont = 15;

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
randompars = make_randomparameters(pars, 0.3);
scalefreepars = make_scalefreeparameters(pars, 3);

f_PRS = figure('Renderer', 'painters', 'Position', [50 800 800 600]);
options = odeset('RelTol', 1.0e-6,'AbsTol', 1.0e-6);

for i = 1:3
    if i == 1
        params = fixeddegreepars;
    elseif i == 2
        params = randompars;
    elseif i == 3
        params = scalefreepars;
    end
    
    % The full scale simulation using the adjacency matrix:
    [t, thetas, A] = DOPRI_simulatenetwork(tnow,tend,IC,h,params);
    z = orderparameter(thetas);
    
    % The OA mean field theory:
    oa_params = prepareOAparameters(params);
    OAIC = ones(oa_params.l,1)*z(1) + 0.001*randn(oa_params.l,1);
    [Toa, b_i] = ode45(@(t,x) MFROA(t,x,oa_params), [tnow, tend], OAIC, options);
    Z_oa = orderparameter_oa(b_i, oa_params.P, oa_params.k, oa_params.N);
    
    % Other order parameters:
    degrees = sum(A,2);
    z_net = orderparameter_net(thetas, degrees, A);
    z_mf = orderparameter_mf(thetas, degrees);
    z_link = orderparameter_link(thetas, degrees, A);

    % Plotting
    subplot(3,1,i); hold on; box on;
    
    plot(t, abs(z), 'LineWidth', 2);
    plot(Toa, abs(Z_oa), 'LineWidth', 2);
    plot(t, abs(z_net), 'LineWidth', 2);
    plot(t, abs(z_mf), 'LineWidth', 2);
    plot(t, abs(z_link), 'LineWidth', 2);

    xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
    ylabel('$\| Z (t) \|$','Interpreter','latex', 'FontSize', labelfont)

    legend('Kuramoto order parameter', 'OA order parameter', 'Network order parameter', 'Mean field order parameter', 'Link field order parameter', 'FontSize', labelfont-5, 'Location', 'southeast')
    removewhitspace();
end







%%




