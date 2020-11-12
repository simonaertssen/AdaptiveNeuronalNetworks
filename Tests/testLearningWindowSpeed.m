clear all; close all; clc
% Test the evaluation speed of some of the learning windows:
%% Setup
N = 1000;
tpts = 30001;

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

%% Observe function shape
dt = linspace(-5.0e-4, 5.0e-4, tpts);

figure; hold on;
plot(dt, Song2000Window(dt));
% plot(dt, Kempter1999Window(dt));
plot(dt, Waddington2014Window(dt))
% plot(dt, ChrolCannon2012Window(dt))
xlabel("Time t [s]"); ylabel("dW")

%% Function properties
clc
integral(@Song2000Window, -5.0e-4, 5.0e-4)
integral(@Waddington2014Window, -5.0e-4, 5.0e-4)
integral(@ChrolCannon2012Window, -100, 100)
integral(@Kempter1999Window, -100, 100)

%% Evaluate multiple time:
clc
tic;
timeme(N, tpts, @Waddington2014Window);
toc;

tic;
timeme(N, tpts, @ChrolCannon2012Window);
toc;

tic;
timeme(N, tpts, @Kempter1999Window);
toc;

% Turns out Waddington is quicker

%% Functions
function out = timeme(N, tpts, handle)
    dt = randn(N,1);
    for i = 1:tpts
        out = handle(dt);
    end
end

function dW = Song2000Window(dt)
    t_pos = 0.020e-3;
    t_neg = 0.020e-3;
    A_pos = 0.005;
    A_neg = 0.00525;
    
    dW = zeros(size(dt));
    t_neg_idx = dt <= 0;
    t_pos_idx = 0 < dt;
    dW(t_neg_idx) = A_pos*exp(dt(t_neg_idx)/t_pos);
    dW(t_pos_idx) = -A_neg*exp(-dt(t_pos_idx)/t_neg);
end

function dW = Waddington2014Window(dt)
    lr = 0.005; alpha = 8.0e-6;
    dW = lr*(1 - ((dt+alpha).^2)./alpha^2).*exp(-abs(dt+alpha)./alpha);
end


function dW = ChrolCannon2012Window(dt)
    Ap = 2.5e-4; 
    Am = 1e-4;
    dW = Ap*exp(-((dt + 15e-7).^2/2e-10)) - Am*exp(-((dt + 15e-7).^2/2e-9));
end


function dW = Kempter1999Window(dt)
    t_syn = 5.0e-6;
    t_pos = 1.0e-6;
    t_neg = 20.0e-6;
    A_p = 1;
    A_n = -1;
    learning_rate = 1.0e-5;
    eps = 1.0e-9;
    
    dW = zeros(size(dt));
    t_neg_idx = dt <= 0;
    t_pos_idx = 0 < dt;
    dW(t_neg_idx) = exp(dt(t_neg_idx)/t_syn + eps).*(A_p*(1-dt(t_neg_idx)/t_pos) + A_n*(1-dt(t_neg_idx)/t_neg));
    dW(t_pos_idx) = A_p*exp(-dt(t_pos_idx)/t_pos + eps) + A_n*exp(-dt(t_pos_idx)/t_neg + eps);
    dW = learning_rate * dW;
end
