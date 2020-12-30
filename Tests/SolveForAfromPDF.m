close all; clear all; clc
addpath('../Functions');

%% Making the adjacency matrix from a degree distribution
% From P(k) ~ k_i and P(k) ~ k_j to A_ij: a combinatorial problem...
% k_in = sum(A,2), k_out = sum(A,1)
% A -> k_in
% |
% v
% k_out
% The challenge is to get k_in(i) for each row, and then to sample from
% k_out so that we get the right number of degrees.
% Try to just always get the largest k values to put in a randomly chosen
% row.

N = 10;
pars.N = N;
networkpars = make_randomparameters(pars, 0.3);

degrees_i = networkpars.degrees_i;
degrees_o = networkpars.degrees_o;

titlefont = 20;
labelfont = 20;

%% Function handle:
f_A = figure('Renderer', 'painters', 'Position', [0 800 1400 400]); 
w = 0.1;
%% Test a fixed degree network:
s = subplot(1, 3, 1); hold on; axis square; box on;
pars.N = 500; vec = linspace(0, pars.N, 11);
netdegree = 100;
fdpars = make_fixeddegreeparameters(pars, netdegree); 

A_fixeddegree = adjacencymatrix(fdpars.degrees_i, fdpars.degrees_o); 

im = imagesc(full(A_fixeddegree)); set(gca,'YDir','reverse');
colormap(gray);
xlim([0, pars.N]); ylim([0, pars.N]);

title('Fixed-degree', 'FontSize', titlefont)
xticks(vec)
xticklabels(string(vec))
% % xlabel('\boldmath$k^{\rm out}$', 'Interpreter', 'latex', 'FontSize', labelfont);
% % ylabel('\boldmath$k^{\rm in}$', 'Interpreter', 'latex', 'FontSize', labelfont);
xlabel('Presynaptic neuron j', 'FontSize', labelfont)
ylabel('Postynaptic neuron i', 'FontSize', labelfont)

% pos = s.Position
% s.Position = [pos(1), pos(2),w, w];

%% Test using the poisson distribution of random networks
subplot(1, 3, 2); hold on; axis square; box on;

netp = 0.200402;
rdpars = make_randomparameters(pars, netp);

A_random = adjacencymatrix(rdpars.degrees_i, rdpars.degrees_o);

im = imagesc(full(A_random)); set(gca,'YDir','reverse');
colormap(gray);
xlim([0, pars.N]); ylim([0, pars.N]);

title('Random', 'FontSize', titlefont)
% xlabel('\boldmath$k^{\rm out}$', 'Interpreter', 'latex', 'FontSize', labelfont);
% ylabel('\boldmath$k^{\rm in}$', 'Interpreter', 'latex', 'FontSize', labelfont);
xlabel('Presynaptic neuron j', 'FontSize', labelfont)
% ylabel('Postynaptic neuron i', 'FontSize', labelfont)

%% Now using scale free networks:
subplot(1, 3, 3); hold on; axis square; box on;

degree = 2.1;
sfpars = make_scalefreeparameters(pars, degree, 50, 260);
sfpars.meandegree

A_scalefree = adjacencymatrix(sfpars.degrees_i, sfpars.degrees_o);

im = imagesc(full(A_scalefree)); set(gca,'YDir','reverse');
colormap(gray)
xlim([0, pars.N]); ylim([0, pars.N]);

title('Scale-free', 'FontSize', titlefont)
% xlabel('\boldmath$k^{\rm out}$', 'Interpreter', 'latex', 'FontSize', labelfont);
% ylabel('\boldmath$k^{\rm in}$', 'Interpreter', 'latex', 'FontSize', labelfont);
xlabel('Presynaptic neuron j', 'FontSize', labelfont)
% ylabel('Postynaptic neuron i', 'FontSize', labelfont)

%% Save:
set(findall(gcf,'-property','FontName'),'FontName','Avenir')
exportgraphics(f_A,'../Figures/Adjacency_matrices.pdf', 'ContentType','vector')
close(f_A)