close all; clear all; clc;
addpath('../Functions');

%% Investigate the properties of different network graphs
pars.N = 100;

if gpuDeviceCount > 0
    d = gpuDevice(gpuDeviceCount);
    disp(d)
end
initarray = make_GPUhandle();

%% Graph properties of a fixed degree network / diracnet:
nsamples = 10;
netdegrees = linspace(log(pars.N), pars.N-1, nsamples); % Use log to make sure k >>> log(N)

Ls = zeros(nsamples,1);
CCs = zeros(nsamples,1);
COs = zeros(nsamples,1);

for i = 1:nsamples
    fdpars = make_fixeddegreeparameters(pars, round(netdegrees(i)));    
    [L, ~, CClosed, ~, COpen, ~] = graphproperties(double(adjacencymatrix(fdpars.degrees_in)));
    Ls(i) = L;
    CCs(i) = CClosed;
    COs(i) = COpen;
    i
end
Ls(Ls == Inf) = 1.0e4;

%% Plotting
netdegreesinterp = linspace(min(netdegrees), max(netdegrees), nsamples^2);

f_fixeddegree = figure('Renderer', 'painters', 'Position', [50 800 600 350]); grid on; hold on;
ylim([0, 5]); title('Fixed degree network', 'FontSize', 20)
plot(netdegreesinterp, interp1(netdegrees,Ls,netdegreesinterp,'pchip'), 'LineWidth', 2)
scatter(netdegrees, Ls, 'xk', 'HandleVisibility','off')

plot(netdegreesinterp, interp1(netdegrees,CCs,netdegreesinterp,'pchip'), 'LineWidth', 2)
scatter(netdegrees, CCs, 'xk', 'HandleVisibility','off')

plot(netdegreesinterp, interp1(netdegrees,COs,netdegreesinterp,'pchip'), 'LineWidth', 2)
scatter(netdegrees, COs, '+k', 'HandleVisibility','off')

xlabel('$$\langle k \rangle$$', 'Interpreter', 'latex', 'FontSize', 15)
legend('$$L(p)$$', '$$C_c(p)$$', '$$C_o(p)$$', 'Location' ,'northeast', 'Interpreter', 'latex', 'FontSize', 20)

removewhitspace();
print(f_fixeddegree, '../Figures/FixedDegreeNetworkProperties.png', '-dpng', '-r300')

%% Graph properties of a random network:
netps = linspace(0.05, 0.95, nsamples); % Use log to make sure k >>> log(N)

Ls = zeros(nsamples,1);
CCs = zeros(nsamples,1);
COs = zeros(nsamples,1);

for i = 1:nsamples
    randompars = make_randomparameters(pars, netps(i));
    [L, ~, CClosed, ~, COpen, ~] = graphproperties(double(adjacencymatrix(randompars.degrees_in)));
    Ls(i) = L;
    CCs(i) = CClosed;
    COs(i) = COpen;
    i
end
Ls(Ls == Inf) = 1.0e2;

%% Plotting
netpsinterp = linspace(min(netps), max(netps), nsamples^2);

f_random = figure('Renderer', 'painters', 'Position', [50 800 600 350]); grid on; hold on;
ylim([0, 5]); title('Fixed degree network', 'FontSize', 20)

plot(netpsinterp, interp1(netps,Ls,netpsinterp,'pchip'), 'LineWidth', 2)
scatter(netps, Ls, 'xk', 'HandleVisibility','off')

plot(netpsinterp, interp1(netps,CCs,netpsinterp,'pchip'), 'LineWidth', 2)
scatter(netps, CCs, 'xk', 'HandleVisibility','off')

plot(netpsinterp, interp1(netps,COs,netpsinterp,'pchip'), 'LineWidth', 2)
scatter(netps, COs, 'xk', 'HandleVisibility','off')

xlabel('$$p$$', 'Interpreter', 'latex', 'FontSize', 15)
legend('$$L(p)$$', '$$C_c(p)$$', '$$C_o(p)$$', 'Location' ,'northeast', 'Interpreter', 'latex', 'FontSize', 20)

removewhitspace();
print(f_random, '../Figures/RandomNetworkProperties.png', '-dpng', '-r300')


%% Graph properties of a random network:
netdegrees = linspace(2.1, 5, nsamples);

Ls = zeros(nsamples,1);
CCs = zeros(nsamples,1);
COs = zeros(nsamples,1);

for i = 1:nsamples
    sfpars = make_scalefreeparameters(pars, netdegrees(i));
    [L, ~, CClosed, ~, COpen, ~] = graphproperties(double(adjacencymatrix(sfpars.degrees_in)));
    Ls(i) = L;
    CCs(i) = CClosed;
    COs(i) = COpen;
    i
end
Ls(Ls == Inf) = 1.0e2;

%% Plotting
netdegreesinterp = linspace(min(netdegrees), max(netdegrees), nsamples^2);

f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 600 350]); grid on; hold on;
ylim([0, 5]); title('Scale free network', 'FontSize', 20)
plot(netdegreesinterp, interp1(netdegrees,Ls,netdegreesinterp,'pchip'), 'LineWidth', 2)
scatter(netdegrees, Ls, 'xk', 'HandleVisibility','off')

plot(netdegreesinterp, interp1(netdegrees,CCs,netdegreesinterp,'pchip'), 'LineWidth', 2)
scatter(netdegrees, CCs, 'xk', 'HandleVisibility','off')

plot(netdegreesinterp, interp1(netdegrees,COs,netdegreesinterp,'pchip'), 'LineWidth', 2)
scatter(netdegrees, COs, '+k', 'HandleVisibility','off')

xlabel('$$\gamma$$', 'Interpreter', 'latex', 'FontSize', 15)
legend('$$L(p)$$', '$$C_c(p)$$', '$$C_o(p)$$', 'Location' ,'northeast', 'Interpreter', 'latex', 'FontSize', 20)

removewhitspace();
print(f_scalefree, '../Figures/ScalefreeNetworkProperties.png', '-dpng', '-r300')
