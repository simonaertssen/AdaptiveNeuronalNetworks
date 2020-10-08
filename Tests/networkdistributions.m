close all; clear all; clc;
addpath('../Functions');

%% Test the different networks and their statistical properties
% We need to evaluate whether some of the random numbers we are pulling are
% actually following the right pdf.

pars.N = 100;

%% A fixed degree network / diracnet: CORRECT
fixeddegreepars = make_fixeddegreeparameters(pars, round(pars.N/2));

limits = round(fixeddegreepars.meandegree + fixeddegreepars.meandegree*[-0.5, 0.5]);

figure; hold on; grid on;
xlim(limits);
x = limits(1):limits(2);
histogram(fixeddegreepars.degrees_in, 'Normalization', 'pdf');
plot(x, fixeddegreepars.P(x)/pars.N);

%% A random network: CORRECT
randompars = make_randomparameters(pars, 0.02);

limits = round(randompars.meandegree + randompars.meandegree*[-1, 1]);

figure; hold on; grid on;
xlim(limits);
x = limits(1):limits(2);
histogram(randompars.degrees_in, 'Normalization', 'pdf');
plot(x, randompars.P(x)/pars.N);

%% A scale free network: CORRECT
kmin = 1000; kmax = 2000;
scalefreepars = make_scalefreeparameters(pars, 3, kmin, kmax);

limits = [kmin, kmax];

figure; hold on; grid on;
x = limits(1):limits(2);
histogram(scalefreepars.degrees_in, 'Normalization', 'pdf');
plot(x, scalefreepars.P(x)/pars.N);


%% A small world network:
A_fixeddegree = adjacencymatrix(fixeddegreepars.degrees_in);

%%
p = 0.3;
idx = find(A_fixeddegree);
freeidx = setdiff(1:pars.N^2, idx);
assert(numel(idx) + numel(freeidx))
disp(numel(idx) + numel(freeidx))

indices = randperm(pars.N, round(0.3*pars.N));






