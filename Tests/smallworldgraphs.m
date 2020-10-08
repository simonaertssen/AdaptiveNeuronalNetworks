close all; clear all; clc;
addpath('../Functions');

%% Make and investigate some small-world networks!
pars.N = 1000;
fixeddegreepars = make_fixeddegreeparameters(pars, round(pars.N/2));

%% A:
A_fixeddegree = adjacencymatrix(fixeddegreepars.degrees_in);

%%
p = 0.3; numrewires = round(0.3*pars.N);
idx = find(A_fixeddegree);

freeidx = setdiff(1:pars.N^2, idx);
assert(numel(idx) + numel(freeidx) == pars.N^2)

idx(randperm(numel(idx), numrewires)) = freeidx(randperm(numel(freeidx), numrewires));

[xidx,yidx] = ind2sub(size(A_fixeddegree),idx);
A = sparse(xidx, yidx, ones(nnz(A_fixeddegree), 1, 'logical'));

immse(full(double(A_fixeddegree)), full(double(A)))

%%
degrees_in = sum(A,2);
degrees_out = sum(A,1); 

histogram(degrees_in, 'Normalization', 'pdf', 'NumBins', round(sqrt(numel(degrees_in))));
