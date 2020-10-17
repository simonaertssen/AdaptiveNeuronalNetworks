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
export = true;

%% Theta model parameters:
tnow = 0; tend = 10;
h = 0.005;

pars.N = 10000;
pars.a_n = 0.666666666666666666667;
seed = 1; rng(seed);
IC = zeros(pars.N, 1); %-pi/2 * ones(pars.N, 1);
IC = -pi/2 * ones(pars.N, 1);

%% Make a GPU init handle:
if gpuDeviceCount > 0
    d = gpuDevice(gpuDeviceCount-1);
    disp(d)
end
initarray = make_GPUhandle();

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
%     A = initarray(adjacencymatrix(params.degrees_in, params.degrees_out));
%     [t, thetas] = ode113(@(t,x,K) thetaneurons_full(t,x,params.K,A,params.e,1/params.meandegree,params.a_n), [tnow, tend], IC, options);
%     thetas = wrapToPi(thetas)';
%     z = orderparameter(thetas);
 
    % The OA mean field theory:
    oa_params = prepareOAparameters(params);
    [TOA, ZOA] = OA_simulatenetwork(tnow, tend, IC, oa_params, options);
    
    % Other order parameters:
    degrees = sum(A,2);
    z_net = orderparameter_net(thetas, degrees, A);
    z_mf = orderparameter_mf(thetas, degrees);
    z_link = orderparameter_link(thetas, degrees, A);

    % Plotting
    imrow(i) = subplot(3,1,i); hold on; box on;
    
    plot(t, abs(z), 'LineWidth', 2);
    plot(TOA, abs(ZOA), 'LineWidth', 2);
    plot(t, abs(z_net), 'LineWidth', 2);
    plot(t, abs(z_mf), 'LineWidth', 2);
    plot(t, abs(z_link), 'LineWidth', 2);
    
    if i < 3
        set(gca, 'XTickLabel', [])
    end
    ylabel('$\| Z (t) \|$','Interpreter','latex', 'FontSize', labelfont)

    removewhitspace();
end

set(imrow(1).Title,'String', sprintf('\\bf Fixed-degree network:  $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fixeddegreepars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex');
set(imrow(2).Title,'String', sprintf('\\bf Random network:  $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, randompars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex');
set(imrow(3).Title,'String', sprintf('\\bf Scale-free network:  $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, scalefreepars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex');

legend('Kuramoto o. p.', 'OA o. p.', 'Network o. p.', 'Mean field o. p.', 'Link field o. p.', 'FontSize', labelfont-5, 'Location', 'southoutside', 'Orientation', 'horizontal')
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
suptitle(sprintf('PSR state: \\eta_0 = %0.1f, \\delta = %0.1f, K = %0.1f', pars.eta0, pars.delta, pars.K))
% exportpdf(f_PRS, '../Figures/InvestigateOrderParametersPRS.pdf', export);
print(f_PRS, '../Figures/InvestigateOrderParametersPRS.pdf', '-dpdf', '-bestfit');
close(f_PRS)

disp('Made PRS state')

%% PSS state: one single stable focus
pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

% Network distributions and parameters:
fixeddegreepars = make_fixeddegreeparameters(pars, round(pars.N*0.3));
randompars = make_randomparameters(pars, 0.3);
scalefreepars = make_scalefreeparameters(pars, 3);

f_PSS = figure('Renderer', 'painters', 'Position', [50 800 800 600]);
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
%     A = initarray(adjacencymatrix(params.degrees_in, params.degrees_out));
%     [t, thetas] = ode113(@(t,x,K) thetaneurons_full(t,x,params.K,A,params.e,1/params.meandegree,params.a_n), [tnow, tend], IC, options);
%     thetas = wrapToPi(thetas)';
%     z = orderparameter(thetas);
 
    % The OA mean field theory:
    oa_params = prepareOAparameters(params);
    [TOA, ZOA] = OA_simulatenetwork(tnow, tend, IC, oa_params, options);
    
    % Other order parameters:
    degrees = sum(A,2);
    z_net = orderparameter_net(thetas, degrees, A);
    z_mf = orderparameter_mf(thetas, degrees);
    z_link = orderparameter_link(thetas, degrees, A);

    % Plotting
    imrow(i) = subplot(3,1,i); hold on; box on;
    
    plot(t, abs(z), 'LineWidth', 2);
    plot(TOA, abs(ZOA), 'LineWidth', 2);
    plot(t, abs(z_net), 'LineWidth', 2);
    plot(t, abs(z_mf), 'LineWidth', 2);
    plot(t, abs(z_link), 'LineWidth', 2);
    
    if i < 3
        set(gca, 'XTickLabel', [])
    end
    ylabel('$\| Z (t) \|$','Interpreter','latex', 'FontSize', labelfont)

    removewhitspace();
end

set(imrow(1).Title,'String', sprintf('\\bf Fixed-degree network:  $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fixeddegreepars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex');
set(imrow(2).Title,'String', sprintf('\\bf Random network:  $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, randompars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex');
set(imrow(3).Title,'String', sprintf('\\bf Scale-free network:  $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, scalefreepars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex');

legend('Kuramoto o. p.', 'OA o. p.', 'Network o. p.', 'Mean field o. p.', 'Link field o. p.', 'FontSize', labelfont-5, 'Location', 'southoutside', 'Orientation', 'horizontal')
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)

suptitle(sprintf('PSS state: \\eta_0 = %0.1f, \\delta = %0.1f, K = %0.1f', pars.eta0, pars.delta, pars.K))
% exportpdf(f_PSS, '../Figures/InvestigateOrderParametersPSS.pdf', export);
print(f_PSS, '../Figures/InvestigateOrderParametersPSS.pdf', '-dpdf', '-bestfit');
close(f_PSS)

disp('Made PSS state')

%% CPW state: the infamous limit cycle
pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

% Network distributions and parameters:
fixeddegreepars = make_fixeddegreeparameters(pars, round(pars.N*0.3));
randompars = make_randomparameters(pars, 0.3);
scalefreepars = make_scalefreeparameters(pars, 3);

f_CPW = figure('Renderer', 'painters', 'Position', [50 800 800 600]);
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
%     A = initarray(adjacencymatrix(params.degrees_in, params.degrees_out));
%     [t, thetas] = ode113(@(t,x,K) thetaneurons_full(t,x,params.K,A,params.e,1/params.meandegree,params.a_n), [tnow, tend], IC, options);
%     thetas = wrapToPi(thetas)';
%     z = orderparameter(thetas);
 
    % The OA mean field theory:
    oa_params = prepareOAparameters(params);
    [TOA, ZOA] = OA_simulatenetwork(tnow, tend, IC, oa_params, options);

    % Other order parameters:
    degrees = sum(A,2);
    z_net = orderparameter_net(thetas, degrees, A);
    z_mf = orderparameter_mf(thetas, degrees);
    z_link = orderparameter_link(thetas, degrees, A);

    % Plotting
    imrow(i) = subplot(3,1,i); hold on; box on;
    
    plot(t, abs(z), 'LineWidth', 2);
    plot(TOA, abs(ZOA), 'LineWidth', 2);
    plot(t, abs(z_net), 'LineWidth', 2);
    plot(t, abs(z_mf), 'LineWidth', 2);
    plot(t, abs(z_link), 'LineWidth', 2);
    
    if i < 3
        set(gca, 'XTickLabel', [])
    end
    ylabel('$\| Z (t) \|$','Interpreter','latex', 'FontSize', labelfont)

    removewhitspace();
end

set(imrow(1).Title,'String', sprintf('\\bf Fixed-degree network:  $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fixeddegreepars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex');
set(imrow(2).Title,'String', sprintf('\\bf Random network:  $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, randompars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex');
set(imrow(3).Title,'String', sprintf('\\bf Scale-free network:  $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, scalefreepars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex');

legend('Kuramoto o. p.', 'OA o. p.', 'Network o. p.', 'Mean field o. p.', 'Link field o. p.', 'FontSize', labelfont-5, 'Location', 'southoutside', 'Orientation', 'horizontal')
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)

suptitle(sprintf('PSS state: \\eta_0 = %0.1f, \\delta = %0.1f, K = %0.1f', pars.eta0, pars.delta, pars.K))
% exportpdf(f_CPW, '../Figures/InvestigateOrderParametersCPW.pdf', export);
print(f_CPW, '../Figures/InvestigateOrderParametersCPW.pdf', '-dpdf', '-bestfit');
close(f_CPW)

disp('Made CPW state')