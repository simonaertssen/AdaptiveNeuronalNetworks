clear all; close all; clc;
% Here we will make phase space plots for the report, illustrating all the
% different global states

%%
addpath('../Functions');
addpath('../Mean Field Reductions/');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')
rect = [-50 50 800 800];
fixeddegreecolor = [0.3010 0.7450 0.9330];
cm = [1,0,0; 0, 0.7410, 0.4470; 0, 0.4470, 0.6410];

titlefont = 15;
labelfont = 50;
linewidth = 4;
arrowst = 12;

fsolveoptions = optimset('Display','off');
odeoptions = odeset('RelTol', 1.0e-8,'AbsTol', 1.0e-8);

export = true;

%% Grids:
% Axis square:
unitsquare = [-1, -1, 1, 1; 1, -1, -1, 1];
th = 0:pi/50:2*pi;
unitcircle = [cos(th); sin(th)];
drawcircle = round(unitcircle);

% Main grid:
l = 1; stp = 2*l/40; interval = -l:stp:l;
[X,Y] = meshgrid(interval,interval);

[in, ~] = inpolygon(X, Y, cos(th), sin(th));
[sz, ~] = size(in); m = ceil(sz/2); l = round(m*0.25);
in(1, :) = 0; in(end, :) = 0; in(:, 1) = 0; in(:, end) = 0;
in(1, m-l:m+l) = 1; in(end, m-l:m+l) = 1; in(m-l:m+l, 1) = 1; in(m-l:m+l, end) = 1;

% Draw grid
% l = 1; stp = 2*l/6; interval = -l:stp:l;
% [X1,Y1] = meshgrid(interval,interval);
% 
% % Startpoints
% [in, on] = inpolygon(X1, Y1, [0.9*cos(th) 0], [0 0.9*sin(th)]);
% [sz, ~] = size(in); m = ceil(sz/2);
% in(1, :) = 0; in(end, :) = 0; in(:, 1) = 0; in(:, end) = 0;
% in(1, m) = 1; in(end, m) = 1; in(m, 1) = 1; in(m, end) = 1;
% X1start = reshape(X1(in), 5, 5); Y1start = reshape(Y1(in), 5, 5);

%% Theta neurons parameters:
pars.N = 5000;
pars.a_n = 0.666666666666666666667;
seed = 2; rng(seed);

%% 1. MFR phase space: PSR
pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

close all;
f_MFRPSR = figure('Renderer', 'painters', 'Position', rect); hold on; box on;

% Equilibrium point:
[eqpt,~,~,~,JACOB] = fsolve(@(x) MFR2D(0,x,pars), [-0.5 -0.5], fsolveoptions);

% Arrows:
complex = MFR(0, X + 1i*Y, pars);
quiver(X, Y, real(complex), imag(complex), 'color', cm(2,:))

% Lines:
complex = MFR(0, X + 1i*Y, pars);
startx = 0.8*cos( -pi/5:pi/5:pi); starty = 0.8*sin(-pi/5:pi/5:pi);
startx(2) = []; starty(2) = [];
scatter(startx, starty, 50, fixeddegreecolor, 'filled', 'o', 'LineWidth', linewidth);
s = streamline(X, Y, real(complex), imag(complex), startx, starty, [0.05,700]);
set(s, 'color', fixeddegreecolor, 'LineWidth', linewidth);
for i = 1:6
    dist = (s(i).XData - eqpt(1)).^2 + (s(i).YData - eqpt(2)).^2 - 0.02;
    dist(end) = -1;
    idx = find(dist < 0, 1, 'first');
    s(i).XData = s(i).XData(1:idx);
    s(i).YData = s(i).YData(1:idx);
    plot_arrow(s(i).XData(end-arrowst), s(i).YData(end-arrowst), s(i).XData(end), s(i).YData(end),'linewidth', 3, ...
    'color', fixeddegreecolor,'facecolor', fixeddegreecolor,'edgecolor', fixeddegreecolor, 'headwidth',2,'headheight',6);
end

phasespaceplot();

% Plot the equilibrium points and Jacobian:
[V,D] = eig(JACOB);
index = 1;
for i = 1:2
    lambda = D(i,i);
    bef = 0.04;
    endl = 0.095;
    draw = 0.75*endl*lambda;
    if lambda < 0
        plot_arrow(eqpt(1) + endl*V(1,index), eqpt(2) + endl*V(2,index), eqpt(1) + bef*V(1,index), eqpt(2) + bef*V(2,index), 'headwidth',0.9,'headheight',6); 
        plot_arrow(eqpt(1) - endl*V(1,index), eqpt(2) - endl*V(2,index), eqpt(1) - bef*V(1,index), eqpt(2) - bef*V(2,index), 'headwidth',0.9,'headheight',6); 
        plot([eqpt(1) - draw*V(1,index), eqpt(1) + draw*V(1,index)], [eqpt(2) - draw*V(2,index), eqpt(2) + draw*V(2,index)], 'k', 'LineWidth', linewidth);
    end
    index = index + 1;
end
scatter(eqpt(1), eqpt(2), 150, 'or', 'filled')

% End figure:
set(gcf,'color','w'); xlim([-1,1]); ylim([-1,1]);
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
xlabel('Re$\left[ Z(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ Z(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
set(findall(gcf,'-property','FontName'),'FontName','Avenir')

exportgraphics(f_MFRPSR,'../Figures/PhaseSpace/MFRPSR.pdf')
close(f_MFRPSR)

%% 2. MFR phase space: PSS
pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

close all
f_MFRPSS = figure('Renderer', 'painters', 'Position', rect); hold on; box on;

% Equilibrium point:
[eqpt,~,~,~,JACOB] = fsolve(@(x) MFR2D(0,x,pars), [-0.5 -0.5], fsolveoptions);

% Arrows:
complex = MFR(0, X + 1i*Y, pars);
quiver(X, Y, real(complex), imag(complex), 'color', cm(2,:))

% Lines:
complex = MFR(0, X + 1i*Y, pars);
startx = 1; starty = 0;
scatter(startx, starty, 50, fixeddegreecolor, 'filled', 'o', 'LineWidth', linewidth);% 'color', fixeddegreecolor);
s = streamline(X, Y, real(complex), imag(complex), startx, starty, [0.05,1550]);
set(s, 'color', fixeddegreecolor, 'LineWidth', linewidth);
plot_arrow(s.XData(end-arrowst), s.YData(end-arrowst), s.XData(end), s.YData(end),'linewidth', 2, ...
    'color', fixeddegreecolor,'facecolor', fixeddegreecolor,'edgecolor', fixeddegreecolor, 'headwidth',0.7,'headheight',3);

phasespaceplot();

% Plot the equilibrium points and Jacobian:
[V,D] = eig(JACOB);
scatter(eqpt(1), eqpt(2), 150, 'or', 'filled')

% End figure:
hold off; set(gcf,'color','w'); xlim([-1,1]); ylim([-1,1]); axis square;
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
xlabel('Re$\left[ Z(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ Z(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

exportgraphics(f_MFRPSS,'../Figures/PhaseSpace/MFRPSS.pdf')
close(f_MFRPSS)

%% 3. MFR phase space: CPW
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);

% Start figure:
close all
f_MFRCPW = figure('Renderer', 'painters', 'Position', rect); hold on; box on;

% Equilibrium point:
[eqpt,~,~,~,~] = fsolve(@(x) MFR2D(0,x,pars), [-0.8 -0.6], fsolveoptions);

% Arrows:
complex = MFR(0, X + 1i*Y, pars);
quiver(X, Y, real(complex), imag(complex), 'color', cm(2,:))

% Lines:
complex = MFR(0, X + 1i*Y, pars);
startx = [0, 0, -0.6, -0.8]; starty = [-1, -0.8, 0.4, 0.2];
scatter(startx, starty, 50, fixeddegreecolor, 'filled', 'o', 'LineWidth', linewidth);
s = streamline(X, Y, real(complex), imag(complex), startx, starty, [0.01,9000]);
set(s, 'color', fixeddegreecolor, 'LineWidth', linewidth);

% Streamlines yield inexact results?
for i = 1:4
    dist = (s(i).XData - eqpt(1)).^2 + (s(i).YData - eqpt(2)).^2 - 0.02;
    dist(end) = -1;
    idx = find(dist < 0, 1, 'first');
    s(i).XData = s(i).XData(1:idx);
    s(i).YData = s(i).YData(1:idx);
    plot_arrow(s(i).XData(end-5*arrowst), s(i).YData(end-5*arrowst), s(i).XData(end), s(i).YData(end),'linewidth', 2, ...
    'color', fixeddegreecolor,'facecolor', fixeddegreecolor,'edgecolor', fixeddegreecolor, 'headwidth',0.7,'headheight',3);     
end

% Limit cycle:
Z0 = -0.2731 - 0.0092*1i;
[T, Z] = ode45(@(t,x) MFR(t,x,pars), [0, 1000], Z0, odeoptions);
Z = flip(Z(round(numel(T)*0.97):end,:));
[~, pksloc] = findpeaks(abs(Z),'MinPeakDistance',100);
idx = pksloc(1):pksloc(3);
plot(real(Z(idx)), imag(Z(idx)), '-', 'LineWidth', linewidth, 'Color', fixeddegreecolor);
plot_arrow(real(Z(idx(end))), imag(Z(idx(end))), real(Z(idx(end)-8)), imag(Z(idx(end)-8)),'linewidth', 2, ...
    'color', fixeddegreecolor,'facecolor', fixeddegreecolor,'edgecolor', fixeddegreecolor, 'headwidth',0.7,'headheight',3);

phasespaceplot();

% Plot the equilibrium points and Jacobian:
ICs = [-0.5 -0.5; -0.8 -0.6];
for eqpointindex = 1:2
    [eqpt,~,~,~,JACOB] = fsolve(@(x) MFR2D(0,x,pars), ICs(eqpointindex,:), fsolveoptions);
    [V,D] = eig(JACOB);
    index = 1;
    for i = 1:2
        lambda = D(i,i);
        bef = 0.04;
        endl = 0.095;
        draw = 0.5*endl*lambda;
        if lambda < 0
            plot_arrow(eqpt(1) + endl*V(1,index), eqpt(2) + endl*V(2,index), eqpt(1) + bef*V(1,index), eqpt(2) + bef*V(2,index), 'headwidth',0.7,'headheight',3);
            plot_arrow(eqpt(1) - endl*V(1,index), eqpt(2) - endl*V(2,index), eqpt(1) - bef*V(1,index), eqpt(2) - bef*V(2,index), 'headwidth',0.7,'headheight',3);
            plot([eqpt(1) - draw*V(1,index), eqpt(1) + draw*V(1,index)], [eqpt(2) - draw*V(2,index), eqpt(2) + draw*V(2,index)], 'k', 'LineWidth', linewidth);
        else
            plot_arrow(eqpt(1) + endl*V(1,index), eqpt(2) + endl*V(2,index),eqpt(1) + draw*V(1,index), eqpt(2) + draw*V(2,index),'headwidth',0.7,'headheight',3);
            plot_arrow(eqpt(1) - endl*V(1,index), eqpt(2) - endl*V(2,index),eqpt(1) - draw*V(1,index), eqpt(2) - draw*V(2,index),'headwidth',0.7,'headheight',3);
            plot([eqpt(1) - draw*V(1,index), eqpt(1) + draw*V(1,index)], [eqpt(2) - draw*V(2,index), eqpt(2) + draw*V(2,index)], 'k', 'LineWidth', linewidth);
        end
        index = index + 1;
    end
    scatter(eqpt(1), eqpt(2), 150, 'or', 'filled')
end

% Draw central focus and draw round vector
[eqpt,FVAL,EXITFLAG,OUTPUT,JACOB] = fsolve(@(x) MFR2D(0,x,pars), [0,0], fsolveoptions);
t = linspace(0, 2*pi + 0.3, 100);
radius = t/max(t)/10;
curve = [eqpt(1) + radius.*cos(t); eqpt(2) + radius.*sin(t)];
plot(curve(1,:), curve(2,:), 'k', 'LineWidth', linewidth)
scatter(eqpt(1), eqpt(2), 150, 'or', 'filled')
plot_arrow(curve(1,end-6), curve(2,end-6), curve(1,end), curve(2,end), 'headwidth',1,'headheight',1.5);

% End figure:
hold off; set(gcf,'color','w'); xlim([-1,1]); ylim([-1,1]); axis square;
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
xlabel('Re$\left[ Z(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ Z(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

exportgraphics(f_MFRCPW,'../Figures/PhaseSpace/MFRCPW.pdf')
close(f_MFRCPW)


%% 4. OA random phase space: PSR
pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_randomparameters(pars, 0.2));

close all
f_OARPSR = figure('Renderer', 'painters', 'Position', rect); hold on; box on;

drawOAvectors(X + 1i*Y, in, p, cm(2,:));

odeoptions = odeset('RelTol', 1.0e-12);
odeoptions.backwards = false;
col = p.colorvec;
startx = 0.8*cos( -pi/5:pi/5:pi); starty = 0.8*sin(-pi/5:pi/5:pi);
startx(2) = []; starty(2) = [];
tlengths = [0.92, 1.15, 1.2, 1.1, 0.85, 0.55];
for i = 1:length(startx)
    OAIC = map_Ztozoa_better(startx(i) + starty(i)*1i, p);
    [~, ZOA] = OA_simulatenetwork(0, tlengths(i), OAIC, p, odeoptions);

    scatter(real(ZOA(1)), imag(ZOA(1)), 50, col, 'filled', 'o', 'LineWidth', linewidth);% 'color', col);

    Zplot = plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth, 'color', col);
    endline = ZOA(end-3) - ZOA(end);
    endpoint = ZOA(end) + 0.04*endline/abs(endline);
    plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
    'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.7,'headheight',3);   
end

% Equilibrium:
eqptb = OA_fixedpointiteration(ones(p.Mk,1), p);
eqptZ = eqptb'*p.P(p.k)/p.N;
scatter(real(eqptZ), imag(eqptZ), 150, 'or', 'filled')

J = MFROAJ(eqptb, p);

phasespaceplot();

% End figure:
hold off; set(gcf,'color','w'); xlim([-1,1]); ylim([-1,1]); axis square;
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

exportgraphics(f_OARPSR,'../Figures/PhaseSpace/MFOARPSR_random.pdf')
close(f_OARPSR)

%% 5. OA random phase space: PSS
pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_randomparameters(pars, 0.3));

close all
f_OARPSS = figure('Renderer', 'painters', 'Position', rect); hold on; box on;
col = p.colorvec;

z0s = drawOAvectors(X + 1i*Y, in, p, cm(2,:));

startx = 1; starty = 0; tlength = 3.4;
odeoptions = odeset('RelTol', 1.0e-12); odeoptions.backwards = false;

OAIC = map_Ztozoa_better(startx + starty*1i, p);
[~, ZOA] = OA_simulatenetwork(0, tlength, OAIC, p, odeoptions);
scatter(real(ZOA(1)), imag(ZOA(1)), 50, col, 'filled', 'o', 'LineWidth', linewidth);% 'color', col);

Zplot = plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth, 'color', col);
endline = ZOA(end-3) - ZOA(end);
endpoint = ZOA(end) + 0.04*endline/abs(endline);
plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.7,'headheight',3); 

phasespaceplot();

% Equilibrium:
eqptb = OA_fixedpointiteration(ones(p.Mk,1), p);
eqptZ = eqptb'*p.P(p.k)/p.N;
scatter(real(eqptZ), imag(eqptZ), 150, 'or', 'filled')

% End figure:
hold off; set(gcf,'color','w'); xlim([-1,1]); ylim([-1,1]); axis square;
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

exportgraphics(f_OARPSS,'../Figures/PhaseSpace/MFOARPSS_random.pdf')
close(f_OARPSS)

%% 6. OA random phase space: CPW
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_randomparameters(pars, 0.3));

close all
f_OARCPW = figure('Renderer', 'painters', 'Position', rect); hold on; box on;
drawfixeddegreelimitcycle();
drawOAvectors(X + 1i*Y, in, p, cm(2,:));

col = p.colorvec;
startx = [0, -0.8, -0.6, 0, 0]; starty = [-0.4, 0.2, 0.4, -1, -0.8];
tlengths = [1.6, 0.5, 0.65, 2.6, 2.15];

for i = 1:length(startx)-1
    OAIC = map_Ztozoa_better(startx(i) + starty(i)*1i, p);
    [~, ZOA] = OA_simulatenetwork(0, tlengths(i), OAIC, p);

    scatter(real(ZOA(1)), imag(ZOA(1)), 50, col, 'filled', 'o', 'LineWidth', linewidth);% 'color', col);

    plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth, 'color', col);
    endline = ZOA(end-3) - ZOA(end);
    endpoint = ZOA(end) + 0.04*endline/abs(endline);
    plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
    'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.7,'headheight',3);   
end

% Limit cycle:
[T, ZOA] = OA_simulatenetwork(0, 100, ones(p.Mk,1)*(-0.73*1i), p);
% plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth, 'color', col);
ZOA = flip(ZOA(round(numel(T)*0.9):end,:));
[~, pksloc] = findpeaks(abs(ZOA),'MinPeakDistance',100);
idx = pksloc(1):pksloc(3);
plot(real(ZOA(idx)), imag(ZOA(idx)), '-', 'LineWidth', linewidth, 'Color', col);
plot_arrow(real(ZOA(end)), imag(ZOA(end)), real(ZOA(end-5)), imag(ZOA(end-5)),'linewidth', 2, ...
    'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.7,'headheight',3);

phasespaceplot();

% Stable quilibrium:
eqptb = OA_fixedpointiteration(ones(p.Mk,1), p);
eqptZ = eqptb'*p.P(p.k)/p.N;
scatter(real(eqptZ), imag(eqptZ), 150, 'or', 'filled')

% End figure:
hold off; set(gcf,'color','w'); xlim([-1,1]); ylim([-1,1]); axis square;
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

exportgraphics(f_OARCPW,'../Figures/PhaseSpace/MFOARCPW_random.pdf')
close(f_OARCPW)


%% 7. OA scalefree phase space: PSR
pars.eta0 = -0.9; pars.delta = 0.8; pars.K = -2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_scalefreeparameters(pars, 3));

close all
f_OARPSR = figure('Renderer', 'painters', 'Position', rect); hold on; box on;

qs = drawOAvectors(X + 1i*Y, in, p, cm(2,:));

odeoptions = odeset('RelTol', 1.0e-12);
odeoptions.backwards = false;
col = p.colorvec;
startx = 0.8*cos( -pi/5:pi/5:pi); starty = 0.8*sin(-pi/5:pi/5:pi);
startx(2) = []; starty(2) = [];
tlengths = [0.92, 1.15, 1.2, 1.1, 0.85, 0.55];

for i = 1:length(startx)
    OAIC = map_Ztozoa_better(startx(i) + starty(i)*1i, p);
    
    [~, ZOA] = OA_simulatenetwork(0, tlengths(i), OAIC, p, odeoptions);

    scatter(real(ZOA(1)), imag(ZOA(1)), 50, col, 'filled', 'o', 'LineWidth', linewidth);% 'color', col);

    Zplot = plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth, 'color', col);
    endline = ZOA(end-3) - ZOA(end);
    endpoint = ZOA(end) + 0.04*endline/abs(endline);
    plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
    'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.7,'headheight',3);   
end
 
% Equilibrium:
eqptb = OA_fixedpointiteration(ones(p.Mk,1), p);
eqptZ = eqptb'*p.P(p.k)/p.N;
scatter(real(eqptZ), imag(eqptZ), 150, 'or', 'filled')

J = MFROAJ(eqptb, p);

phasespaceplot();

% End figure:
hold off; set(gcf,'color','w'); xlim([-1,1]); ylim([-1,1]); axis square;
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
exportgraphics(f_OARPSR,'../Figures/PhaseSpace/MFOARPSR_scalefree.pdf')
close(f_OARPSR)

%% 8. OA scalefree phase space: PSS
pars.eta0 = 0.5; pars.delta = 0.7; pars.K = 2;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_scalefreeparameters(pars, 3, 750, 1000));

close all
f_OARPSS = figure('Renderer', 'painters', 'Position', rect); hold on; box on;
col = p.colorvec;

z0s = drawOAvectors(X + 1i*Y, in, p, cm(2,:));

startx = 1; starty = 0; tlength = 3.4;
odeoptions = odeset('RelTol', 1.0e-12); odeoptions.backwards = true;

[~, ZOA] = OA_simulatenetwork(0, tlength, ones(p.Mk,1)*(startx + starty*1i), p, odeoptions);
scatter(startx, starty, 50, col, 'filled', 'o', 'LineWidth', linewidth);% 'color', col);
Zplot = plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth, 'color', col);
endline = ZOA(end-3) - ZOA(end);
endpoint = ZOA(end) + 0.04*endline/abs(endline);
plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.7,'headheight',3); 

phasespaceplot();

% Equilibrium:
eqptb = OA_fixedpointiteration(ones(p.Mk,1), p);
eqptZ = eqptb'*p.P(p.k)/p.N;
scatter(real(eqptZ), imag(eqptZ), 150, 'or', 'filled')

% End figure:
hold off; set(gcf,'color','w'); xlim([-1,1]); ylim([-1,1]); axis square;
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

exportgraphics(f_OARPSS,'../Figures/PhaseSpace/MFOARPSS_scalefree.pdf')
close(f_OARPSS)

%% 9. OA scalefree phase space: CPW
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
p = prepareOAparameters(make_scalefreeparameters(pars, 3));

close all
f_OARCPW = figure('Renderer', 'painters', 'Position', rect); hold on; box on;
drawfixeddegreelimitcycle();
z0s = drawOAvectors(X + 1i*Y, in, p, cm(2,:));

odeoptions = odeset('RelTol', 1.0e-6); odeoptions.backwards = false;
col = p.colorvec;

startx = [-0.8, 0, 0]; starty = [0.4, -1, -0.8];
tlengths = [0.5, 0.65, 3, 1];
bw = -0.5;
for i = 1:length(startx)
%     OAIC = map_Ztozoa_better(startx(i) + starty(i)*1i, p);
    OAIC = ones(p.Mk,1)*(startx(i) + starty(i)*1i);
    
    [~, ZOA] = OA_simulatenetwork(0, tlengths(i), OAIC, p, odeoptions);
    
    scatter(real(ZOA(1)), imag(ZOA(1)), 50, col, 'filled', 'o', 'LineWidth', linewidth);

    Zplot = plot(real(ZOA), imag(ZOA), 'LineWidth', linewidth, 'color', col);
    endline = ZOA(end-3) - ZOA(end);
    endpoint = ZOA(end) + 0.04*endline/abs(endline);
    plot_arrow(real(endpoint), imag(endpoint), real(ZOA(end)), imag(ZOA(end)),'linewidth', 2, ...
    'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.7,'headheight',3);   
end

% Limit cycle:
[T, ZOA] = OA_simulatenetwork(0, 100, ones(p.Mk,1)*(-0.2*1i), p, odeoptions);
ZOA = flip(ZOA(round(numel(T)*0.9):end,:));
[~, pksloc] = findpeaks(abs(ZOA),'MinPeakDistance',100);
idx = pksloc(1):pksloc(3);
plot(real(ZOA(idx)), imag(ZOA(idx)), '-', 'LineWidth', linewidth, 'Color', col);
plot_arrow(real(ZOA(end)), imag(ZOA(end)), real(ZOA(end-3)), imag(ZOA(end-3)),'linewidth', 2, ...
    'color', col,'facecolor', col,'edgecolor', col, 'headwidth',0.7,'headheight',3);

phasespaceplot();

% Stable quilibrium:
eqptb = OA_fixedpointiteration(ones(p.Mk,1), p);
eqptZ = eqptb'*p.P(p.k)/p.N;
scatter(real(eqptZ), imag(eqptZ), 150, 'or', 'filled')

% End figure:
hold off; set(gcf,'color','w'); xlim([-1,1]); ylim([-1,1]); axis square;
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

exportgraphics(f_OARCPW,'../Figures/PhaseSpace/MFOARCPW_scalefree.pdf')
close(f_OARCPW)

%% 10. Difference in the limit cycles of the scale-free network
tnow = 0; tend = 50;
h = 0.01;

pars.N = 5000;
pars.a_n = 0.666666666666666666667;
pars.eta0 = 10.75; pars.delta = 0.5; pars.K = -9;

seed = 2; rng(seed); IC = wrapToPi(randn(pars.N, 1)*1.4);
pars.e = randcauchy(seed, pars.eta0, pars.delta, pars.N);
odeoptions = odeset('RelTol', 1.0e-9,'AbsTol', 1.0e-9);

degree = 3;
IC = pi*ones(pars.N, 1) - pi;

sfpars = make_scalefreeparameters(pars, degree);
[~, thetasfull, A] = DOPRI_simulatenetwork(tnow,tend,IC,h,sfpars);
zfull = orderparameter(thetasfull);
ts = findlimitcycle(abs(zfull));
zfull = zfull(ts(1):ts(2));
disp('Full scale test done')

sfpars = prepareOAparameters(sfpars);
z0 = orderparameter(IC)*ones(sfpars.Mk,1);
[~, ZOA] = OA_simulatenetwork(tnow, tend, z0, sfpars, odeoptions);
TOAs = findlimitcycle(abs(ZOA));
ZOA = ZOA(TOAs(1):TOAs(2));
disp('OA mean field test done')

% Figure:
labelfont = 30;
col = sfpars.colorvec;
f_scalefree = figure('Renderer', 'painters', 'Position', rect); hold on; box on;
cycle = drawfixeddegreelimitcycle();
cycle.HandleVisibility = 'off';
z0s = drawOAvectors(X + 1i*Y, in, sfpars, cm(2,:));

plot(real(zfull), imag(zfull), '-', 'LineWidth', 3, 'Color', '#0072BD');
plot(real(ZOA), imag(ZOA), '-', 'LineWidth', 3, 'Color', col);
    
phasespaceplot();

legend('$$Z(t)_{A_{ij}}$$', '$$\bar{Z}(t)_{MF_{OA}}$$', 'Interpreter', 'latex', 'FontSize', labelfont, 'Location', 'northeast', 'Orientation','horizontal')

% End figure:
hold off; set(gcf,'color','w'); xlim([-1,1]); ylim([-1,1]); axis square;
set(gca,'XTickLabel',[]); set(gca,'YTickLabel',[]); 
xlabel('Re$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)
ylabel('Im$\left[ \bar{Z}(t)\right]$','Interpreter','latex', 'FontSize', labelfont)

exportgraphics(f_scalefree,'../Figures/PhaseSpace/ScalefreeLimCycles.pdf')
close(f_scalefree)
