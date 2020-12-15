% In this script we will be studying the reponse of the theta neuron model
% on different types of input current.

clear all; close all; clc;

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

tnow = 0; tend = 10;
F = @thetaneuron; h = 0.001;

fI = figure('Renderer', 'painters', 'Position', [50 800 1000 200]);
titlefont = 15;
labelfont = 20;
axesfont = 17;
m = 1; n = 3;

%% Investigate the frequency - current curve:
% We know that T = pi/sqrt(I)

imrow(1) = subplot(m,n,1); hold on; box on;
I = @fI_current;
h = 0.01;

tend = 50;

% Theoratical result:
Idraw = linspace(-10,10, 200);
plot(Idraw, sqrt(Idraw)./pi, 'LineWidth', 6, 'color', '#EDB120');

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

%% Draw the explanation about the PRC:
I = 5;
t = linspace(0, 2, 1001);
realthetas = @(t) 2*atan(-sqrt(I)*cot(t*sqrt(I)));
drawrealthetas = spikesNaN(realthetas(t));
realspiketime = t(isnan(drawrealthetas));

phasebifftime = 0.8; change = 0.2;
bifthetas = realthetas(t);
bifidx = find(t == phasebifftime);
bifthetas(bifidx:end) = realthetas(t(bifidx:end) - change);
drawbifthetas = spikesNaN(bifthetas);
bifspiketime = t(isnan(drawbifthetas));

replacedvalue = bifthetas(t == phasebifftime);
bifthetas(bifidx-1) = NaN;
drawbifthetas(bifidx-1) = NaN;

newvaluetime = fsolve(@(t) realthetas(t) - replacedvalue, phasebifftime, optimoptions('fsolve','Display','off'));

imrow(2) = subplot(m,n,2); hold on; box on;
hold on
line([phasebifftime, phasebifftime], [drawrealthetas(bifidx), drawbifthetas(bifidx)], 'LineStyle', '--', 'color', 'k', 'LineWidth', 1)

plot(t, bifthetas, ':k', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, drawbifthetas, '-r', 'LineWidth', 2);

plot(t, realthetas(t), ':k', 'LineWidth', 1, 'HandleVisibility','off');
plot(t, drawrealthetas, '-', 'LineWidth', 2, 'color', '#0072BD');

text(phasebifftime + 0.02, 0.5*(drawbifthetas(bifidx+1) + drawrealthetas(bifidx-1)) + 0.4, '$$\varepsilon$$', 'HorizontalAlignment', 'left', 'FontSize', labelfont, 'Interpreter', 'latex');


xlabel('Time','Interpreter','tex','FontSize', labelfont, 'FontName', 'SansSerif')
ylabel('$\theta$','Interpreter','latex', 'FontSize', labelfont)
set(gca,'XTick',[0, newvaluetime, phasebifftime, realspiketime, bifspiketime], 'XTickLabel',{'0','$$\phi_{\rm{new}}$$','$$\phi$$','$$T$$','$$T_{\phi}$$'}, 'TickLabelInterpreter', 'latex', 'FontSize', axesfont, 'FontName', 'SansSerif')
set(gca,'YTick',[-pi, -pi/2, replacedvalue, 0, pi/2, pi], 'YTickLabel',{'-$$\pi$$','$$-\frac{\pi}{2}$$',' ','0','$$\frac{\pi}{2}$$','$$\pi$$'}, 'TickLabelInterpreter', 'latex', 'YLim', [-pi-0.2, pi + 0.2], 'FontSize', axesfont)


%% Draw the PRC itself:
imrow(3) = subplot(m,n,3); hold on; box on;

I = 2; T = pi/sqrt(I);
t = linspace(0, T, 100);
t = t(1:end-1);
iPRC = 1/(2*I) * (1 - cos(2*pi*t/T));
realthetas = @(t) 2*atan(-sqrt(I)*cot(t*sqrt(I)));
realdthetas = @(t) 2*I*csc(t*sqrt(I)).^2./(I*cot(t*sqrt(I)).^2 + 1);

xlim([0, T]);
yyaxis right;
plot(t, iPRC, 'LineWidth', 2, 'color', '#D95319')
ylabel('\sl PRC', 'FontSize', labelfont, 'FontName', 'SansSerif')

yyaxis left;
plot(t, realthetas(t), 'LineWidth', 2, 'color', '#0072BD')
ylabel('$$\theta$$', 'FontSize', labelfont, 'Interpreter', 'latex')
set(gca,'YTick',[-pi, -pi/2, 0, pi/2, pi], 'YTickLabel',{'-$$\pi$$','$$-\frac{\pi}{2}$$','0','$$\frac{\pi}{2}$$','$$\pi$$'}, 'TickLabelInterpreter', 'latex', 'YLim', [-pi-0.2, pi + 0.2], 'FontSize', axesfont)
ax = gca; ax.YAxis(1).Color = 'k'; ax.YAxis(2).Color = '#D95319';


% es = logspace(-2, -1, 5)
% es = linspace(0.001, 0.05, 5);
% phi = t;
% 
% for e = es
%     Tf = 1/sqrt(I)*(pi/2 - atan(e/sqrt(I) + cot(phi*sqrt(I)))) - phi;
%         PRC = - Tf*(2*I);
%     plot(phi, PRC, 'LineWidth', 2)
% end

xlabel('$\phi$','Interpreter','latex', 'FontSize', labelfont)

set(gca,'XTick',[0, T/2, T], 'XTickLabel',{'0','$$T/2$$','$$T$$'}, 'TickLabelInterpreter', 'latex')


%%
exportgraphics(fI,'../Figures/ThetaNeuronfIandPRC.pdf')


%% Functions:

function I = fI_current(t)
    I = 1.0e-3.*t.^(1);
end