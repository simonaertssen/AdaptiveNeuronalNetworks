% We will test whether the analytical solution 
% theta = 2*np.arctan(np.sqrt(i)*np.tan(t*np.sqrt(i) - np.pi/2))
% is valid by recording the error between simullations and theory.

clear all; close all; clc;

%% Setup
addpath('../Functions');

set(groot,'DefaultAxesXGrid','on')
set(groot,'DefaultAxesYGrid','on')

tnow = 0; tend = 10;
F = @thetaneuron; h = 0.001;

titlefont = 15;
labelfont = 15;

%% 
fexcite = figure('Renderer', 'painters', 'Position', [50 800 1000 200]);


theta = 2*np.arctan(np.sqrt(i)*np.tan(t*np.sqrt(i) - np.pi/2))

%% Functions
function I = simplecurrent(t)
    I = ones(size(t));
end

function I = parabolacurrent(t)
    I = t*t;
end

function I = sinecurrent(t)
    I = sin(t);
end