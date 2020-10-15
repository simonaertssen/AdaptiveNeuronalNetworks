close all; clear all; clc;
addpath('../Functions');

%% Make and investigate some small-world networks!
% Whe rewiring, the in-degree has to change! Outgoing connections stay the
% same.
pars.N = 1000;
fixeddegreepars = make_fixeddegreeparameters(pars, 10);

%% Get fixed degree matrix:
A_fixeddegree = adjacencymatrix(fixeddegreepars.degrees_in);

%% Interchange a percentage of links
A_fixeddegree = adjacencymatrix(fixeddegreepars.degrees_in);

A = A_fixeddegree;
idx = find(A);
[xidx,yidx] = ind2sub(size(A),idx);

p = 0.12; 
rewires = rand(numel(idx), 1) < p;
numrewires = sum(rewires);
rewireidx = find(rewires);

notallowed_idx = [idx; sub2ind(size(A), 1:pars.N, 1:pars.N)'];

freeidx = setdiff(1:pars.N^2, notallowed_idx);
assert(numel(freeidx) + numel(notallowed_idx) == pars.N^2)
[freexidx,freeyidx] = ind2sub(size(A),freeidx);

for j = 1:numrewires
    yidx(rewireidx(j)) = randsample(freeyidx(freeyidx == yidx(rewireidx(j))), 1);
end

% delete rewired links:
A(idx) = 0;
idx = sub2ind(xidx,yidx);
A(idx) = 1;

diffcols = full(sum(A_fixeddegree,2)) - full(sum(A,2));
diffrows = full(sum(A_fixeddegree,1)) - full(sum(A,1));

sum(diffcols)
sum(diffrows)
N2 = pars.N^2; thresh = 1.0e-9;


%%
degrees_in = sum(A,2);
degrees_out = sum(A,1); 

histogram(degrees_in, 'Normalization', 'pdf', 'NumBins', round(sqrt(numel(degrees_in))));

%% Strogatz1998
% Do we obtain the same behaviour as in Strogatz1998? Let's see.
pars.N = 500;
fixeddegreepars = make_fixeddegreeparameters(pars, 10);
A_fixeddegree = adjacencymatrix(fixeddegreepars.degrees_in);

nsamples = 10;
ps = logspace(-4, 0, nsamples);
Ls = zeros(nsamples,1);
CCs = zeros(nsamples,1);
COs = zeros(nsamples,1);

[L0, ~, CC0, ~, CO0, ~] = graphproperties(double(A_fixeddegree));

for i = 1:nsamples
    p = ps(i); 
    A = A_fixeddegree;
    idx = find(A);
    [xidx,yidx] = ind2sub(size(A),idx);

    rewires = rand(numel(idx), 1) < p;
    numrewires = sum(rewires);
    rewireidx = find(rewires);

    notallowed_idx = [idx; sub2ind(size(A), 1:pars.N, 1:pars.N)'];

    freeidx = setdiff(1:pars.N^2, notallowed_idx);
    assert(numel(freeidx) + numel(notallowed_idx) == pars.N^2)
    [freexidx,freeyidx] = ind2sub(size(A),freeidx);

    for j = 1:numrewires
        yidx(rewireidx(j)) = randsample(freeyidx(freeyidx == yidx(rewireidx(j))), 1);
    end

    % delete rewired links:
    A(idx(rewires)) = 0;
    idx = sub2ind(xidx,yidx);
    A(idx) = 1;

    [L, ~, CClosed, ~, COpen, ~] = graphproperties(double(A));
    Ls(i) = L;
    CCs(i) = CClosed;
    COs(i) = COpen;
    i
end

%%
pinterp = logspace(-4, 0, nsamples^2);
figure; grid on; hold on
plot(pinterp, interp1(ps,Ls/L0,pinterp), 'LineWidth', 2)
plot(pinterp, interp1(ps,CCs/CC0,pinterp), 'LineWidth', 2)
%plot(pinterp, interp1(ps,COs/CO0,pinterp), 'LineWidth', 2)
set(gca, 'XScale', 'log')
xlabel('$$p$$', 'Interpreter', 'latex')
legend('$$L(p)/L_0$$', '$$C(p)/C_0$$', 'Location' ,'southwest', 'Interpreter', 'latex', 'FontSize', 20)
% Saved as 'nosmallworldfromdirac.png'


%% Old:
p = 0.1; numrewires = sum(rand(pars.N));
idx = find(A_fixeddegree);

freeidx = setdiff(1:pars.N^2, idx);
assert(numel(idx) + numel(freeidx) == pars.N^2)

% Link a certain number and delink another part: 
rewireidx = randperm(numel(idx), numrewires*2);
idx(rewireidx(1:numrewires)) = freeidx(randperm(numel(freeidx), numrewires));
idx(rewireidx((numrewires+1):end)) = 0;

[xidx,yidx] = ind2sub(size(A_fixeddegree),idx);
A = sparse(xidx, yidx, ones(nnz(A_fixeddegree), 1, 'logical'));

immse(full(double(A_fixeddegree)), full(double(A)))

%% Check
degrees_in = sum(A,2);
degrees_out = sum(A,1); 

histogram(degrees_in, 'Normalization', 'pdf', 'NumBins', round(sqrt(numel(degrees_in))));

%% Strogatz1998
% Do we obtain the same behaviour as in Strogatz1998? Let's see.
pars.N = 500;
fixeddegreepars = make_fixeddegreeparameters(pars, 100);
A_fixeddegree = adjacencymatrix(fixeddegreepars.degrees_in, fixeddegreepars.degrees_out);

nsamples = 10;
ps = logspace(-4, 0, nsamples);
Ls = zeros(nsamples,1);
CCs = zeros(nsamples,1);
COs = zeros(nsamples,1);

[L0, ~, CC0, ~, CO0, ~] = graphproperties(double(A_fixeddegree));

for i = 1:nsamples
    p = ps(i); 
    numrewires = round(p*pars.N);
    
    idx = find(A_fixeddegree);
    freeidx = setdiff(1:pars.N^2, idx);
    assert(numel(idx) + numel(freeidx) == pars.N^2)

    % Link a certain number and delink another part: 
    rewireidx = randperm(numel(idx), numrewires*2);
    idx(rewireidx(1:numrewires)) = freeidx(randperm(numel(freeidx), numrewires));

    [xidx,yidx] = ind2sub(size(A_fixeddegree),idx);
    A = sparse(xidx, yidx, ones(nnz(A_fixeddegree), 1));
    A(idx(rewireidx((numrewires+1):end))) = 0;

    [L, ~, CClosed, ~, COpen, ~] = graphproperties(A);
    Ls(i) = L;
    CCs(i) = CClosed;
    COs(i) = COpen;
    i
end

%%
pinterp = logspace(-4, 0, nsamples^2);
f = figure; grid on; hold on
plot(pinterp, interp1(ps,Ls/L0,pinterp), 'LineWidth', 2)
plot(pinterp, interp1(ps,CCs/CC0,pinterp), 'LineWidth', 2)
%plot(pinterp, interp1(ps,COs/CO0,pinterp), 'LineWidth', 2)
set(gca, 'XScale', 'log')
xlabel('$$p$$', 'Interpreter', 'latex')
legend('$$L(p)/L_0$$', '$$C(p)/C_0$$', 'Location' ,'southwest', 'Interpreter', 'latex', 'FontSize', 20)
% Saved as 'nosmallworldfromdirac.png'
exportpdf(export, f, 'nosmallworldfromdirac');

%% Realisation:
% We don't need to change elements in the adjacency matrix to model the
% displaced links, we can just 'add' a vector of zero-mean.
% If p*N links will change, p*N links will be added and p*N links will be
% removed.

linkstoadd = randfixedsum(pars.N, 1, numrewires, 0, numrewires);
sum(linkstoadd)
