clear all; close all; clc
% Test the speed of some of the learning windows and measure some more
% properties. Make graphs.

%% Setup
addpath('../Functions');

N = 1000;
tpts = 30001;

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

titlefont = 15;
labelfont = 13;
export = true;

%% Observe function shapes: biphasic windows
dt = linspace(-0.08, 0.1, tpts);

f_windows = figure('Renderer', 'painters', 'Position', [50 800 600 200]); 

subplot(1,2,1); box on; hold on; xlim([dt(1), dt(end)]);
plot(dt, Kempter1999Window(dt), 'LineWidth', 2, 'Color', '#0072BD')
legend("$$W(t)_K$$", 'Interpreter', 'latex', 'FontSize', labelfont, 'Orientation','horizontal', 'Location', 'southeast')
xlabel("$t$ [s]", 'Interpreter', 'latex', 'FontSize', labelfont); 
ylabel("$W$", 'Interpreter', 'latex', 'FontSize', labelfont);

subplot(1,2,2); ylim([-0.14, 0.12]); xlim([dt(1), dt(end)]); box on; hold on;
yleft = Song2000Window(dt(dt<=0)); yright = Song2000Window(dt(dt>0));
plot(dt(dt<=0), yleft, 'LineWidth', 2, 'Color', '#D95319'); 
plot(dt(dt>0), yright, 'LineWidth', 2, 'Color', '#D95319', 'HandleVisibility', 'off')
line([0, 0],[min(yleft), max(yright)],'Color','k','LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off')
scatter(0, Song2000Window(0), 30, [0.8500 0.3250 0.0980], 'filled', 'o')

legend("$$W(t)_S$$", 'Interpreter', 'latex', 'FontSize', labelfont, 'Orientation','horizontal', 'Location', 'southeast')
xlabel("$t$ [s]", 'Interpreter', 'latex', 'FontSize', labelfont); 
ylabel("$W$", 'Interpreter', 'latex', 'FontSize', labelfont);

% exportpdf(f_windows, '../Figures/Learning/LearningWindowsBiphasic.pdf', export);
print(f_windows, '../Figures/Learning/LearningWindowsBiphasic.png', '-dpng', '-r400')

close(f_windows)

%% Observe function shapes: triphasic windows
dt = linspace(-0.08, 0.1, tpts);

f_windows = figure('Renderer', 'painters', 'Position', [50 800 600 200]); box on; hold on;
subplot(1,2,1); box on; hold on; xlim([dt(1), dt(end)]); ylim([-0.14, 0.12]);
y = ChrolCannon2012Window(dt);
plot(dt, y, 'LineWidth', 2, 'Color', '#77AC30')
legend("$$W(t)_C$$", 'Interpreter', 'latex', 'FontSize', labelfont, 'Orientation','horizontal', 'Location', 'southwest')
xlabel("$t$ [s]", 'Interpreter', 'latex', 'FontSize', labelfont); 
ylabel("$W$", 'Interpreter', 'latex', 'FontSize', labelfont);

subplot(1,2,2); box on; hold on; xlim([dt(1), dt(end)]); ylim([-0.14, 0.12]);
plot(dt, Waddington2014Window(dt), 'LineWidth', 2, 'Color', '#A2142F')

legend("$$W(t)_W$$", 'Interpreter', 'latex', 'FontSize', labelfont, 'Orientation','horizontal', 'Location', 'southwest')
xlabel("$t$ [s]", 'Interpreter', 'latex', 'FontSize', labelfont); 
ylabel("$W$", 'Interpreter', 'latex', 'FontSize', labelfont);

% exportpdf(f_windows, '../Figures/Learning/LearningWindowsTriphasic.pdf', export);
print(f_windows, '../Figures/Learning/LearningWindowsTriphasic.png', '-dpng', '-r400')
close(f_windows)


%% Observe function shapes: three first windows
dt = linspace(-0.08, 0.1, tpts);

f_windows = figure('Renderer', 'painters', 'Position', [50 800 1000 200]); 

subplot(1,3,1); box on; hold on; xlim([dt(1), dt(end)]);
plot(dt, Kempter1999Window(dt), 'LineWidth', 2, 'Color', '#0072BD')
legend("$$W(t)_K$$", 'Interpreter', 'latex', 'FontSize', labelfont, 'Orientation','horizontal', 'Location', 'southeast')
xlabel("$t$ [s]", 'Interpreter', 'latex', 'FontSize', labelfont); 
ylabel("$W$", 'Interpreter', 'latex', 'FontSize', labelfont);

subplot(1,3,2); ylim([-0.14, 0.12]); xlim([dt(1), dt(end)]); box on; hold on;
yleft = Song2000Window(dt(dt<=0)); yright = Song2000Window(dt(dt>0));
plot(dt(dt<=0), yleft, 'LineWidth', 2, 'Color', '#D95319'); 
plot(dt(dt>0), yright, 'LineWidth', 2, 'Color', '#D95319', 'HandleVisibility', 'off')
line([0, 0],[min(yleft), max(yright)],'Color','k','LineStyle','--', 'LineWidth', 1, 'HandleVisibility', 'off')
scatter(0, Song2000Window(0), 30, [0.8500 0.3250 0.0980], 'filled', 'o')

legend("$$W(t)_S$$", 'Interpreter', 'latex', 'FontSize', labelfont, 'Orientation','horizontal', 'Location', 'southeast')
xlabel("$t$ [s]", 'Interpreter', 'latex', 'FontSize', labelfont); 
% ylabel("$W$", 'Interpreter', 'latex', 'FontSize', labelfont);

subplot(1,3,3); box on; hold on; xlim([dt(1), dt(end)]); ylim([-0.14, 0.12]);
y = ChrolCannon2012Window(dt);
plot(dt, y, 'LineWidth', 2, 'Color', '#77AC30')
legend("$$W(t)_C$$", 'Interpreter', 'latex', 'FontSize', labelfont, 'Orientation','horizontal', 'Location', 'southwest')
xlabel("$t$ [s]", 'Interpreter', 'latex', 'FontSize', labelfont); 
% ylabel("$W$", 'Interpreter', 'latex', 'FontSize', labelfont);

% print(f_windows, '../Figures/Learning/LearningWindows.png', '-dpng', '-r400')
exportgraphics(f_windows,'../Figures/Learning/LearningWindows.pdf')

close(f_windows)

%% Function properties
clc
integral(@Kempter1999Window, -0.08, 0.1)
integral(@Song2000Window, -0.08, 0.1)
integral(@Song2017Window, -0.08, 0.1)
integral(@ChrolCannon2012Window, -0.08, 0.1)
integral(@Waddington2014Window, -0.08, 0.1)

%% Evaluate multiple time:
clc

tic;
timeme(N, 1000, @Kempter1999Window);
toc;

tic;
timeme(N, 1000, @Song2000Window);
toc;

tic;
timeme(N, 1000, @ChrolCannon2012Window);
toc;

tic;
timeme(N, 1000, @Waddington2014Window);
toc;

% Turns out Waddington is quicker


%% Observe function shapes: ip learning
dt = linspace(0, 0.25, tpts);

f_windows = figure('Renderer', 'painters', 'Position', [50 800 300 200]); box on; hold on;
y = Song2017ip(dt);
xline(0.09, '--k')
xline(0.11, '--k')

yleft = Song2017ip(dt(dt<0.09)); ymiddle = Song2017ip(dt((dt>0.09) & (dt<0.11))); yright = Song2017ip(dt(dt>0.11));
plot(dt(dt<0.09), yleft, 'LineWidth', 2, 'Color', '#0072BD'); 
plot(dt((dt>0.09) & (dt<0.11)), ymiddle, 'LineWidth', 2, 'Color', '#0072BD', 'HandleVisibility', 'off'); 
plot(dt(dt>0.11), yright, 'LineWidth', 2, 'Color', '#0072BD', 'HandleVisibility', 'off')
scatter(0.09, Song2017ip(0.09), 30, [0 0.4470 0.7410], 'filled', 'o')
scatter(0.11, Song2017ip(0.11), 30, [0 0.4470 0.7410], 'filled', 'o')

text(0.089, -0.04, '$$T_{ \rm min}$$', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right')
text(0.112, -0.04, '$$T_{ \rm max}$$', 'Interpreter', 'latex', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left')

xlabel("$t$ [s]", 'Interpreter', 'latex', 'FontSize', labelfont); 
ylabel("$\phi_i(t)$", 'Interpreter', 'latex', 'FontSize', labelfont);

print(f_windows, '../Figures/IPlearningFunction.png', '-dpng', '-r400')
close(f_windows)


%% Functions
function out = timeme(N, tpts, handle)
    dt = randn(N,1);
    for i = 1:tpts
        out = handle(dt);
    end
end

function dW = Kempter1999Window(dt)
    t_syn = 5;
    t_pos = 1;
    t_neg = 20;
    A_p = 1;
    A_n = -1;
    learning_rate = 1.0e-5;
    eps = 1.0e-9;
    dt = -dt * 1.0e3; % Convert to seconds
    
    dW = zeros(size(dt));
    t_neg_idx = dt <= 0;
    t_pos_idx = 0 < dt;
    dW(t_neg_idx) = exp(dt(t_neg_idx)/t_syn + eps).*(A_p*(1-dt(t_neg_idx)/t_pos) + A_n*(1-dt(t_neg_idx)/t_neg));
    dW(t_pos_idx) = A_p*exp(-dt(t_pos_idx)/t_pos + eps) + A_n*exp(-dt(t_pos_idx)/t_neg + eps);
    dW = learning_rate * dW;
end

function dW = Kempter1999WindowOld(dt)
    t_syn = 5;
    t_pos = 1;
    t_neg = 20;
    A_p = 1;
    A_n = -1;
    learning_rate = 5.0e-2;
    eps = 1.0e-9;
    dt = -dt * 1.0e3; % Convert to seconds
    
    dW = zeros(size(dt));
    t_neg_idx = dt <= 0;
    t_pos_idx = 0 < dt;
    dW(t_neg_idx) = exp(dt(t_neg_idx)/t_syn + eps).*(A_p*(1-dt(t_neg_idx)/t_pos) + A_n*(1-dt(t_neg_idx)/t_neg));
    dW(t_pos_idx) = A_p*exp(-dt(t_pos_idx)/t_pos + eps) + A_n*exp(-dt(t_pos_idx)/t_neg + eps);
    dW = learning_rate * dW;
end


function dW = Song2000Window(dt)
    t_pos = 20;
    t_neg = 20;
    A_pos = 0.1;
    A_neg = -0.12;

    dt = dt * 1.0e3; % Convert to seconds
    dW = zeros(size(dt));
    t_zero = dt == 0;
    t_neg_idx = dt <= 0;
    t_pos_idx = 0 < dt;
    dW(t_zero) = 0;
    dW(t_pos_idx) = A_pos*exp(-dt(t_pos_idx)/t_pos);
    dW(t_neg_idx) = A_neg*exp(dt(t_neg_idx)/t_neg);
end


function dW = Song2017Window(dt)
    t_pos = 20;
    t_neg = 20;
    A_pos = 0.005;
    A_neg = -0.00525;

    dt = dt * 1.0e3; % Convert to seconds
    dW = zeros(size(dt));
    t_pos_idx = dt > 0;
    t_neg_idx = dt <= 0;
    dW(t_pos_idx) = A_pos*exp(-dt(t_pos_idx)/t_pos);
    dW(t_neg_idx) = A_neg*exp(dt(t_neg_idx)/t_neg);
end


function dW = Waddington2014Window(dt)
%     lr = 0.005; alpha = 8.0e-6;
    lr = 0.1; alpha = 4.0;
    dt = dt * 1.0e3; % Convert to seconds
    dW = lr*(1 - ((dt-alpha).^2)./alpha^2).*exp(-abs(dt-alpha)./alpha);
end


function dW = ChrolCannon2012Window(dt)
%     Ap = 2.3e-4; 
%     Am = 1.5e-4;
%     dW = Ap*exp(-((dt + 15e-7).^2/2e-10)) - Am*exp(-((dt + 20e-7).^2/2e-9));

    Ap = 0.23; 
    Am = 0.15;
    t_pos = 200;
    t_neg = 2000;
    dt = dt * 1.0e3; % Convert to seconds
    dW = Ap*exp(-((dt - 15).^2/t_pos)) - Am*exp(-((dt - 20).^2/t_neg));
end


function phi = Song2017ip(dt)
    T_min = 90; %ms
    T_max = 110; %ms
    lr = 0.012;

    dt = dt * 1.0e3; % Convert to seconds
    phi = zeros(size(dt));
    t_neg_idx = dt < T_min;    
    t_pos_idx = dt > T_max;
    
    phi(t_neg_idx) = -lr*exp((T_min - dt(t_neg_idx))/T_min);
    phi(t_pos_idx) = lr*exp((dt(t_pos_idx) - T_max)/T_max);
    end
