close all; clear vars; clc;
addpath('../Functions');

%% Test the different networks and their statistical properties
% We need to evaluate whether some of the random numbers we are pulling are
% actually following the right pdf.

pars.N = 10000;

%% A fixed degree network / diracnet: CORRECT
fixeddegreepars = make_fixeddegreeparameters(pars, 100);

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

%% A scale free network:
kmin = 1000; kmax = 2000;
scalefreepars = make_scalefreeparameters(pars, 3, kmin, kmax);

limits = [kmin, kmax];

figure; hold on; grid on;
x = limits(1):limits(2);
histogram(scalefreepars.degrees_in, 'Normalization', 'pdf');
plot(x, scalefreepars.P(x)/pars.N);


%% A scale free network:

%% Setup: the pdf as a benchmark
addpath('../Functions');
rng(0);

N = 2000; mu = -5; gamma = 2.7;

xloc = mu + linspace(-10*gamma, 10*gamma, N);
pdffunc = @(x) 1 ./ (pi*gamma*(1 + power((x - mu)/gamma, 2)));
pdf = pdffunc(xloc);

%% 1. Using matlab 'makedist':
tic;
pd = makedist('tLocationScale','mu',mu,'sigma',gamma,'nu',1);
e1 = random(pd, N, 1);
toc
% Elapsed time is 0.003010 seconds.

%% 2. Using 'cauchy' on Matlab File Exchange
tic;
r = rand(N,1);
r(r < 0 | 1 < r) = NaN;
e2 = mu + gamma.*tan(pi*(r-0.5));
e2(r == 0)=	-Inf;
e2(r == 1)=	Inf;
toc
% Elapsed time is 0.002249 seconds.

%% 3. Using simple inverse:
tic;
x0 = 1.e-6;
e3 = reshape(mu + gamma*tan(pi*(linspace(x0, 1-x0, N) - 0.5)), N, 1);
toc
% Elapsed time is 0.001055 seconds.

%% Plotting:
f = figure; hold on
plot(xloc, pdf, 'LineWidth', 2.5);
xlim([mu - 10*gamma, mu + 10*gamma])

h1 = histogram(e1, 'BinEdges', xloc(1:50:end), 'Normalization', 'pdf');
h2 = histogram(e2, 'BinEdges', xloc(1:50:end), 'Normalization', 'pdf');
h3 = histogram(e3, 'BinEdges', xloc(1:50:end), 'Normalization', 'pdf');

h1values = h1.Values;
h2values = h2.Values;
h3values = h3.Values;

% close(f)

%% Error:
geterror = @(pdf, hist) sum(power(pdf(1:50:end-50) - hist, 2));
geterror(pdf, h1values) % = 0.0015
geterror(pdf, h2values) % = 0.0013
geterror(pdf, h3values) % = 2.7637e-04



