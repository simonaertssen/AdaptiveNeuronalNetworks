% In this script we will be studying the reponse of the theta neuron model
% on different types of input current.

clear all; close all; clc;

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

tnow = 0; tend = 10;
F = @thetaneuron; h = 0.001;

fexcite = figure('Renderer', 'painters', 'Position', [50 800 1000 200]);
titlefont = 15;
labelfont = 20;
axesfont = 15;
m = 1; n = 3;

%% Excitability:
% The theta neuron model is of type 1 excitability, meaning that the frequency 
% of spikes goes up as the input current increases.
I = @excitabilitycurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
NaNthetas = spikesNaN(thetas);

imrow(1) = subplot(m,n,1); hold on; box on;

%title("Class 1 excitability", 'FontSize', titlefont, 'FontName', 'SansSerif');
yyaxis left; hold on;
plot(t, thetas, ':k', 'LineWidth', 1);
plot(t, NaNthetas, '-', 'LineWidth', 2, 'color', '#0072BD');
ylabel('$\theta$','Interpreter','latex', 'FontSize', labelfont);
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
set(gca,'YTick',-pi:pi:pi, 'YTickLabel',{'-\pi','0','\pi'}, 'TickLabelInterpreter', 'tex', 'YLim', [-pi - 1.5, pi + 0.2], 'FontSize', axesfont, 'FontName','Avenir')
    
yyaxis right; hold on; box on;
maxy = max(I(t));
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca,'YTick', 0:maxy:maxy, 'YLim', [-10, maxy*8], 'FontName','Avenir');


%% Spiking behaviour:
% The neuron spikes when the input current > 0
I = @spikecurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
NaNthetas = spikesNaN(thetas);

imrow(2) = subplot(m,n,2); hold on; box on;

%title("Spiking", 'FontSize', titlefont, 'FontName', 'SansSerif');
yyaxis left; hold on;
plot(t, thetas, ':k', 'LineWidth', 1);
plot(t, NaNthetas, '-', 'LineWidth', 2, 'color', '#0072BD');
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
set(gca,'YTick',-pi:pi:pi, 'YTickLabel',{'-\pi','0','\pi'}, 'TickLabelInterpreter', 'tex', 'YLim', [-pi - 1.5, pi + 0.2], 'FontSize', axesfont, 'FontName','Avenir')

yyaxis right; hold on; box on;
maxy = max(I(t));
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca,'YTick', 0:maxy:maxy, 'YLim', [-2, maxy*8]);


%% Bursting behaviour:
% The neuron spikes when the input current > 0
I = @burstcurrent;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, -pi, h, I);
NaNthetas = spikesNaN(thetas);

imrow(3) = subplot(m,n,3); hold on; box on;

%title("Bursting", 'FontSize', titlefont, 'FontName', 'SansSerif');
yyaxis left; hold on;
plot(t, thetas, ':k', 'LineWidth', 1);
plot(t, NaNthetas, '-', 'LineWidth', 2, 'color', '#0072BD');
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)
%set(gca,'YTick',-pi:pi/2:pi, 'YTickLabel',{'-$$\pi$$','-$$\frac{\pi}{2}$$','0','$$\frac{\pi}{2}$$','$$\pi$$'}, 'TickLabelInterpreter', 'latex', 'YLim', [-pi - 1.5, pi + 0.2], 'FontSize', axesfont)
set(gca,'YTick',-pi:pi:pi, 'YTickLabel',{'-\pi','0','\pi'}, 'TickLabelInterpreter', 'tex', 'YLim', [-pi - 1.5, pi + 0.2], 'FontSize', axesfont, 'FontName','Avenir')

yyaxis right; hold on; box on;
maxy = max(I(t));
plot(t, I(t), '-', 'LineWidth', 1, 'color', '#A2142F');
ylim([-1, max(I(t))*5]);
ylabel('$I$','Interpreter','latex', 'FontSize', labelfont)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#A2142F';
set(gca,'YTick', 0:maxy:maxy, 'YLim', [-200, maxy*8])


%% Save the figure:
set(findall(gcf,'-property','FontName'),'FontName','Avenir')
exportgraphics(fexcite,'../Figures/ThetaNeuronResponseToCurrent.pdf')

%% Investigate the frequency - current curve:
% We know that T = pi/sqrt(I)
% 
fI = figure('Renderer', 'painters', 'Position', [50 800 500 200]); hold on; box on;
I = @fI_current;
h = 0.01;

tend = 50;

% Theoratical result:
Idraw = linspace(-10,10, 200);
plot(Idraw, sqrt(Idraw)./pi, 'LineWidth', 6);

nmeasure = 33;
Imeasure = [linspace(-10,-1, nmeasure/3), linspace(0,1, nmeasure/3), linspace(2,10, nmeasure/3)];
frequencies = zeros(nmeasure,1);
for i = 1:nmeasure
    I = @(t) Imeasure(i);

    % Simulate
    [t, thetas] = DOPRI_singleneuron(F, 0, tend, -pi, h, I);
    NaNthetas = spikesNaN(thetas);

    % Measurements:
    if sum(isnan(NaNthetas)) > 0
        idx = [1, find(isnan(NaNthetas))];
        frequencies(i) = mean(1./diff(t(idx))) + 0;
    else    
        frequencies(i) = 0;
    end
end

drawfrequencies = interp1(Imeasure',frequencies,Idraw);

plot(Imeasure, frequencies, ':k', 'LineWidth', 3);

xlabel('$I$','Interpreter','latex', 'FontSize', labelfont)
ylabel('$\frac{1}{T}$','Interpreter','latex','Rotation', 0, 'FontSize', labelfont*1.5, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
legend('$\frac{\sqrt{I}}{\pi}$', '$1/\hat{T}$', 'Interpreter','latex', 'FontSize', labelfont, 'Location', 'northwest')

% exportgraphics(fI,'../Figures/ThetaNeuronResponseToCurrentPeriod.pdf')

%% Functions:
function I = excitabilitycurrent(t)
    I = max(t/10, power(t,2));
end

function I = spikecurrent(t)
    I = zeros(size(t));
    spike = 10;
    I(t > 1 & t < 2) = spike;
    I(t > 3 & t < 4) = spike;
    
    I(t > 5 & t < 6) = spike;
    I(t > 7 & t < 8) = spike;
end

function I = burstcurrent(t)
    I = zeros(size(t));
    burst = 500;
    I(t > 2 & t < 3) = burst;
    I(t > 6 & t < 7) = burst;
end

function I = fI_current(t)
    I = 1.0e-3.*t.^(1);
end
