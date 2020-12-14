% Display the different analytical solutions
clear all; close all; clc;

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

tnow = 0; tend = 20;
F = @thetaneuron; h = 0.001; IC = -pi;

%%
feqpts = figure('Renderer', 'painters', 'Position', [50 800 1000 200]); 
titlefont = 15;
labelfont = 15;
axesfont = 15;
m = 1; n = 3;


%% Excitable regime
I = @(x) -5;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, IC, h, I);
drawthetas = spikesNaN(thetas);

sI = sqrt(-I(0));
realthetas = 2*atan(2*sI./(1 - exp(2*t*sI)) - sI);

imrow(1) = subplot(m,n,1); hold on; box on;

ylim([-pi - 1.0, pi + 0.1]);
plot(t, thetas, ':k', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, drawthetas, '-', 'LineWidth', 2, 'color', '#0072BD');

plot(t, realthetas, ':k', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, spikesNaN(realthetas), '-', 'LineWidth', 2, 'color', '#77AC30');

ylabel('$\theta$','Interpreter','latex', 'FontSize', labelfont);
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)

legend('Simulation', 'Analytical');


%% Bifurcation
I = @(x) 0;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, IC, h, I);
drawthetas = spikesNaN(thetas);

sI = sqrt(-I(0));
realthetas = 2*atan(2*sI./(1 - exp(2*t*sI)) - sI);

imrow(2) = subplot(m,n,2); hold on; box on;

ylim([-pi - 1.0, pi + 0.1]);
plot(t, thetas, ':k', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, drawthetas, '-', 'LineWidth', 2, 'color', '#0072BD');

plot(t, realthetas, ':k', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, spikesNaN(realthetas), '-', 'LineWidth', 2, 'color', '#77AC30');

ylabel('$\theta$','Interpreter','latex', 'FontSize', labelfont);
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)

legend('Simulation', 'Analytical');

%% Periodic regime
I = @(x) 5;

[t, thetas] = DOPRI_singleneuron(F, tnow, tend, IC, h, I);
drawthetas = spikesNaN(thetas);

sI = sqrt(-I(0));
realthetas = 2*atan(2*sI./(1 - exp(2*t*sI)) - sI);

imrow(3) = subplot(m,n,3); hold on; box on;

ylim([-pi - 1.0, pi + 0.1]);
plot(t, thetas, ':k', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, drawthetas, '-', 'LineWidth', 2, 'color', '#0072BD');

plot(t, realthetas, ':k', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, spikesNaN(realthetas), '-', 'LineWidth', 2, 'color', '#77AC30');

ylabel('$\theta$','Interpreter','latex', 'FontSize', labelfont);
xlabel('$t$','Interpreter','latex', 'FontSize', labelfont)

legend('Simulation', 'Analytical');


%% Save the figure:
print(fanalytic, '../Figures/ThetaNeuronAnalyticalSolution.png', '-dpng', '-r300')


%% Functions

function I = simplecurrent(t)
    I = 2*ones(size(t));
end

function I = linearcurrent(e, t)
    I = 2 + e.*t;
end

function I = paraboliccurrent(t)
    I = 0.05.*t.*t;
end

function I = sinecurrent(t)
    I = 1.1 + sin(t./100);
end