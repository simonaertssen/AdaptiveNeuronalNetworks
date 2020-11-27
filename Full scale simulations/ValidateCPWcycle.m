clear all; close all; clc;
% Simulate a full scale fixed degree network and test whether the results
% are correct, with respect to the order parameters.

%% Setup
addpath('../Functions');
addpath('../Mean Field Reductions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 15;
export = false;

%% Make a GPU init handle:
if gpuDeviceCount > 0
    d = gpuDevice(3);
end
initarray = make_GPUhandle();

%% Theta model parameters:
tnow = 0; tend = 20;
h = 0.005;

pars.N = 10000;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;

seed = 2; rng(seed);
IC = wrapToPi(randn(pars.N, 1)*1.4);

pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
odeoptions = odeset('RelTol', 1.0e-9,'AbsTol', 1.0e-9);

% %% 0. Perform a full scale simulation of a FULLY CONNECTED network:
% % The full scale simulation using the adjacency matrix:
% fdpars = make_fixeddegreeparameters(pars, pars.N);
% [~, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,fdpars);
% zfull = orderparameter(thetasfull);
% ts = findlimitcycle(abs(zfull));
% disp('Full scale test done')
% 
% % The mean field theory for fixed degree networks:
% [~, Z] = ode45(@(t,x) MFR2(t,x,pars), [tnow, tend], gather(zfull(1)), odeoptions);
% Ts = findlimitcycle(abs(Z));
% disp('Mean field test done')
% 
% % The OA mean field theory:
% fdpars = prepareOAparameters(fdpars);
% z0 = map_thetatozoa(gather(thetasfull(:,1)), fdpars);
% [~, ZOA] = OA_simulatenetwork(tnow, tend, gather(z0), fdpars, odeoptions);
% TOAs = findlimitcycle(abs(ZOA));
% disp('OA mean field test done')
% 
% zfull = zfull(ts(1):ts(2));
% Z = Z(Ts(1):Ts(2));
% ZOA = ZOA(TOAs(1):TOAs(2));
% 
% %% Plotting the results:
% f_fullyconnected = figure('Renderer', 'painters', 'Position', [50 800 500 500]); box on; hold on; axis square;
% cycle = drawfixeddegreelimitcycle();
% cycle.HandleVisibility = 'off';
% 
% scatter(real(zfull(1)), imag(zfull(1)), 50, 'x', 'MarkerEdgeColor', '#0072BD', 'LineWidth', 1, 'HandleVisibility', 'off');
% plot(real(zfull), imag(zfull), '-', 'LineWidth', 2, 'Color', '#0072BD');
% 
% scatter(real(Z(1)), imag(Z(1)), 50, '+', 'MarkerEdgeColor', '#D95319', 'LineWidth', 1, 'HandleVisibility', 'off');
% plot(real(Z), imag(Z), '-', 'LineWidth', 2, 'Color', '#D95319');
% 
% scatter(real(ZOA(1)), imag(ZOA(1)), 50, 'o', 'MarkerEdgeColor', '#000000', 'LineWidth', 1, 'HandleVisibility', 'off');
% plot(real(ZOA), imag(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');
% 
% xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
% ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
% 
% phasespaceplot();
% 
% title(sprintf('\\bf Fully connected network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')
% legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southoutside', 'Orientation','horizontal')
% exportpdf(f_fullyconnected, '../Figures/InspectMeanFieldFullyConnectedPhaseSpace.pdf', export);
% close(f_fullyconnected)
% 
% disp('Made fully connected network figure')
% 
% %% 1. Perform a full scale simulation of a fixed degree network:
% netdegree = round(pars.N*0.5);
% 
% % The full scale simulation using the adjacency matrix:
% fdpars = make_fixeddegreeparameters(pars, netdegree);
% [~, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,fdpars);
% zfull = orderparameter(thetasfull);
% ts = findlimitcycle(abs(zfull));
% disp('Full scale test done')
% 
% % The mean field theory for fixed degree networks:
% [~, Z] = ode45(@(t,x) MFR(t,x,pars), [tnow, tend], gather(zfull(1)), odeoptions);
% Ts = findlimitcycle(abs(Z));
% disp('Mean field test done')
% 
% % The OA mean field theory:
% fdpars = prepareOAparameters(fdpars);
% z0 = map_thetatozoa(gather(thetasfull(:,1)), fdpars);
% [~, ZOA] = OA_simulatenetwork(tnow, tend, z0, fdpars, odeoptions);
% TOAs = findlimitcycle(abs(ZOA));
% disp('OA mean field test done')
% 
% 
% zfull = zfull(ts(1):ts(2));
% Z = Z(Ts(1):Ts(2));
% ZOA = ZOA(TOAs(1):TOAs(2));
% 
% %% Plotting the results:
% f_fixeddegree = figure('Renderer', 'painters', 'Position', [50 800 500 500]); box on; hold on; axis square;
% 
% cycle = drawfixeddegreelimitcycle();
% cycle.HandleVisibility = 'off';
% 
% scatter(real(zfull(1)), imag(zfull(1)), 50, 'x', 'MarkerEdgeColor', '#0072BD', 'LineWidth', 1, 'HandleVisibility', 'off');
% plot(real(zfull), imag(zfull), '-', 'LineWidth', 2, 'Color', '#0072BD');
% 
% scatter(real(Z(1)), imag(Z(1)), 50, '+', 'MarkerEdgeColor', '#D95319', 'LineWidth', 1, 'HandleVisibility', 'off');
% plot(real(Z), imag(Z), '-', 'LineWidth', 2, 'Color', '#D95319');
% 
% scatter(real(ZOA(1)), imag(ZOA(1)), 50, 'o', 'MarkerEdgeColor', '#000000', 'LineWidth', 1, 'HandleVisibility', 'off');
% plot(real(ZOA), imag(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');
% 
% xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
% ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
% 
% phasespaceplot();
% 
% title(sprintf('\\bf Fixed degree network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, fdpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')
% legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southoutside', 'Orientation','horizontal')
% exportpdf(f_fixeddegree, '../Figures/InspectMeanFieldFixedDegreePhaseSpace.pdf', export);
% close(f_fixeddegree)
% 
% disp('Made fixed degree network figure')
% 
% %% 2. Perform a full scale simulation of a random network:
% % The full scale simulation using the adjacency matrix:
% netp = 0.3;
% rdpars = make_randomparameters(pars, netp);
% [~, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,rdpars);
% zfull = orderparameter(thetasfull);
% ts = findlimitcycle(abs(zfull));
% disp('Full scale test done')
% 
% % The OA mean field theory:
% rdpars = prepareOAparameters(rdpars);
% z0 = map_thetatozoa(gather(thetasfull(:,1)), rdpars);
% [~, ZOA] = OA_simulatenetwork(tnow, tend, z0, rdpars, odeoptions);
% TOAs = findlimitcycle(abs(ZOA));
% disp('OA mean field test done')
% 
% 
% zfull = zfull(ts(1):ts(2));
% ZOA = ZOA(TOAs(1):TOAs(2));
% 
% %% Plotting the results:
% f_random = figure('Renderer', 'painters', 'Position', [50 800 500 500]); box on; hold on; axis square;
% 
% cycle = drawfixeddegreelimitcycle();
% cycle.HandleVisibility = 'off';
% 
% scatter(real(zfull(1)), imag(zfull(1)), 50, 'x', 'MarkerEdgeColor', '#0072BD', 'LineWidth', 1, 'HandleVisibility', 'off');
% plot(real(zfull), imag(zfull), '-', 'LineWidth', 2, 'Color', '#0072BD');
% 
% scatter(real(ZOA(1)), imag(ZOA(1)), 50, 'o', 'MarkerEdgeColor', '#000000', 'LineWidth', 1, 'HandleVisibility', 'off');
% plot(real(ZOA), imag(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');
% 
% xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
% ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
% 
% phasespaceplot();
% 
% title(sprintf('\\bf Random network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$p$$ = %0.1f', pars.N, rdpars.meandegree, rdpars.netp), 'FontSize', titlefont, 'Interpreter', 'latex')
% legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southoutside', 'Orientation','horizontal')
% exportpdf(f_random, '../Figures/InspectMeanFieldRandomPhaseSpace.pdf', export);
% close(f_random)
% 
% disp('Made random network figure')
% 
%% 3. Perform a full scale simulation of a scale-free network:
% The full scale simulation using the adjacency matrix:
degree = 3;
IC = wrapToPi(randn(pars.N, 1)*2.0);
IC = pi*ones(pars.N, 1) - pi;

sfpars = make_scalefreeparameters(pars, degree);

%%
vec = linspace(sfpars.kmin, sfpars.kmax, sfpars.kmax-sfpars.kmin+1);
[x,y] = meshgrid(vec, vec);
idx = 1:numel(vec);

sfpars.P = @(x,y) P2D(x, y, sfpars.degree, pars.N);

Xin = sfpars.kmin:10:sfpars.kmax;
Yin = sfpars.kmin:10:sfpars.kmax;

Xmat = ones(length(Yin),1)*Xin;
Ymat = Yin'*ones(1,length(Xin));

for i = 1:pars.N
    [sfpars.degrees_i(i),sfpars.degrees_o(i)] = pinky(Xin, Yin, P2D(Xmat, Ymat, degree, pars.N));
end
% sfpars.degrees_o = sfpars.degrees_i(randperm(pars.N));


% Remove extra connections:
% n = sum(sfpars.degrees_i) - sum(sfpars.degrees_o);
% if n > 0
%     abs(
% else
% end

% figure; hold on;
% x = sfpars.kmin:sfpars.kmax;
% plot(x, sfpars.P(x)/sfpars.N)
% histogram(sfpars.degrees_i, 'Normalization', 'pdf')
% histogram(sfpars.degrees_o, 'Normalization', 'pdf')

%%
% sfpars = prepareOAparameters2DP(sfpars);
% 
% figure; hold on;
% surf(x(idx,idx),y(idx,idx),sfpars.P(x(idx,idx),y(idx,idx))/pars.N);
% kminkmax = linspace(sfpars.kmin, sfpars.kmax, 10);
% histogram2(sfpars.degrees_i, sfpars.degrees_o, kminkmax, kminkmax, 'Normalization', 'pdf'); % Normal degree vectors from before


%%
[~, thetasfull, A] = DOPRI_simulatenetwork(tnow,tend,IC,h,sfpars);
zfull = orderparameter(thetasfull);
% ts = findlimitcycle(abs(zfull));
disp('Full scale test done')

% The OA mean field theory:
sfpars = prepareOAparameters2DP(sfpars);
sfpars.Mk
pars.N
% z0 = map_thetatozoa(gather(thetasfull(:,1)), sfpars);
z0 = orderparameter(IC)*ones(sfpars.Mk,1);
[~, ZOA] = OA_simulatenetwork(tnow, tend, z0, sfpars, odeoptions);
% TOAs = findlimitcycle(abs(ZOA));
disp('OA mean field test done')


% zfull = zfull(ts(1):ts(2));
% ZOA = ZOA(TOAs(1):TOAs(2));

%% Plotting the results:
f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 500 500]); box on; hold on; axis square;

cycle = drawfixeddegreelimitcycle();
cycle.HandleVisibility = 'off';

scatter(real(zfull(1)), imag(zfull(1)), 50, 'x', 'MarkerEdgeColor', '#0072BD', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(real(zfull), imag(zfull), '-', 'LineWidth', 2, 'Color', '#0072BD');

scatter(real(ZOA(1)), imag(ZOA(1)), 50, 'o', 'MarkerEdgeColor', '#000000', 'LineWidth', 1, 'HandleVisibility', 'off');
plot(real(ZOA), imag(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');

xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

phasespaceplot();

title(sprintf('\\bf Scale-free network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f, $$\\gamma$$ = %0.1f', pars.N, sfpars.meandegree, sfpars.degree), 'FontSize', titlefont, 'Interpreter', 'latex')
legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southoutside', 'Orientation','horizontal')
% exportpdf(f_scalefree, '../Figures/InspectMeanFieldScaleFreePhaseSpace.pdf', export);
print(f_scalefree, '../Figures/testScaleFree.png', '-dpng', '-r300')
% close(f_scalefree)

disp('Made scale-free network figure')

% %% 4. Perform a full scale simulation of a lognorm network:
% % The full scale simulation using the adjacency matrix:
% lnpars = make_lognormparameters(pars, 3, 1, round(pars.N/5));
% [~, thetasfull] = DOPRI_simulatenetwork(tnow,tend,IC,h,lnpars);
% zfull = orderparameter(thetasfull);
% ts = findlimitcycle(abs(zfull));
% disp('Full scale test done')
% 
% % The OA mean field theory:
% lnpars = prepareOAparameters(lnpars);
% z0 = map_thetatozoa(gather(thetasfull(:,1)), lnpars);
% [~, ZOA] = OA_simulatenetwork(tnow, tend, z0, lnpars, odeoptions);
% TOAs = findlimitcycle(abs(ZOA));
% disp('OA mean field test done')
% 
% zfull = zfull(ts(1):ts(2));
% ZOA = ZOA(TOAs(1):TOAs(2));
% 
% %% Plotting the results:
% f_lognorm = figure('Renderer', 'painters', 'Position', [50 800 500 500]); box on; hold on; axis square;
% 
% cycle = drawfixeddegreelimitcycle();
% cycle.HandleVisibility = 'off';
% 
% scatter(real(zfull(1)), imag(zfull(1)), 50, 'x', 'MarkerEdgeColor', '#0072BD', 'LineWidth', 1, 'HandleVisibility', 'off');
% plot(real(zfull), imag(zfull), '-', 'LineWidth', 2, 'Color', '#0072BD');
% 
% scatter(real(ZOA(1)), imag(ZOA(1)), 50, 'o', 'MarkerEdgeColor', '#000000', 'LineWidth', 1, 'HandleVisibility', 'off');
% plot(real(ZOA), imag(ZOA), '-', 'LineWidth', 2, 'Color', '#000000');
% 
% xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
% ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
% 
% phasespaceplot();
% 
% title(sprintf('\\bf Lognorm network: $$N$$ = %d, $$\\langle k \\rangle$$ = %0.1f', pars.N, lnpars.meandegree), 'FontSize', titlefont, 'Interpreter', 'latex')
% legend('$$Z(t)_{A_{ij}}$$', '$$\overline{Z(t)}_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'southoutside', 'Orientation','horizontal')
% exportpdf(f_lognorm, '../Figures/InspectMeanFieldLogNormPhaseSpace.pdf', export);
% close(f_lognorm)
% 
% disp('Made lognorm network figure')
% 

function P = P2D(x,y,degree,N)
    P = x.^(-degree) + y.^(-degree);
    P = P/sum(P, 'all')*N;
end


function p = prepareOAparameters2DP(p)
    [p.k, ~, ic] = unique([p.degrees_i, p.degrees_o], 'rows', 'stable');
    p.kcount = accumarray(ic, 1);
    p.Mk = numel(p.k)/2;
    
    p.OA = zeros(p.Mk, p.Mk);
    for i = 1:p.Mk
        p.OA(i, :) = p.P(p.k(:,1), p.k(:,2)).*assortativity(p.k(:,1), p.k(:,2), p.k(i,1), p.k(i,2), p.N, p.meandegree, 0);
    end
    p.OA = p.K*p.OA/p.meandegree;
end
