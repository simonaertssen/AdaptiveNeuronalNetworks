clear all; close all; clc
% Test the evaluation speed of some of the learning windows:
%% Setup
N = 10000;
tpts = 30000;

%% Observe function shape
dt = linspace(-100, 100, tpts);
figure; hold on;
plot(dt, Waddington2014Window(dt))
plot(dt, ChrolCannon2012Window(dt))

%% Function properties
integral(@Waddington2014Window, -100, 100)
integral(@ChrolCannon2012Window, -100, 100)

%% Evaluate multiple time:
clc
tic;
timeme(N, tpts, @Waddington2014Window);
toc;

tic;
timeme(N, tpts, @ChrolCannon2012Window);
toc;

% Turns out Waddington is quicker

%% Functions
function out = timeme(N, tpts, handle)
    dt = randn(N,1);
    for i = 1:tpts
        out = handle(dt);
    end
end

function dW = Waddington2014Window(dt)
    lr = 0.1; alpha = 4.0;
    dW = lr*(1 - ((dt-alpha).^2)./alpha^2).*exp(-abs(dt-alpha)./alpha);
end


function dW = ChrolCannon2012Window(dt)
    Ap = 0.25; Am = 0.1;
    dW = Ap*exp(-((dt - 15).^2/200)) - Am*exp(-((dt - 15).^2/2000));
end