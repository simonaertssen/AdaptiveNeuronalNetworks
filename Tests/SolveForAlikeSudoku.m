close all; clear vars; clc
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

N = 1000;

%% Test a fixed degree network:
netdegree = 600;
degrees_in = zeros(N,1);
degrees_in(randperm(N)) = netdegree;
degrees_out = degrees_in(randperm(N));

assert(sum(degrees_in) == sum(degrees_out))

A_fixeddegree = adjacencymatrix(degrees_in, degrees_out);
f_fixeddegree = figure('Renderer', 'painters', 'Position', [50 800 400 400]);
imshow(full(A_fixeddegree));
print(f_fixeddegree, '../Figures/A_fixeddegree.png', '-dpng', '-r300')

%% Test using the poisson distribution of random networks
netp = 0.2;
meandegree = netp*(N - 1);
numlinks = netp*N*(N-1)/2;

degrees_in = poissrnd(meandegree,N,1);
degrees_out = degrees_in(randperm(N)); % because these two have to contain the same number of elements

assert(sum(degrees_in) == sum(degrees_out))

A_random = adjacencymatrix(degrees_in, degrees_out);
f_random = figure('Renderer', 'painters', 'Position', [50 800 400 400]);
imshow(full(A_random));
print(f_random, '../Figures/A_random.png', '-dpng', '-r300')

%% Now using scale free networks:
degree = 3;
kmin = 750; kmax = 2000;

P = @(x) scalefreepdf(x, N, degree, kmin, kmax);

x = linspace(kmin,kmax,N);

degrees_in = randsample(kmin + P(x)', N, true);
degrees_out = degrees_in(randperm(N));
histogram(degrees_in, 'Normalization', 'pdf')
%%
assert(sum(degrees_in) == sum(degrees_out))

A_scalefree = adjacencymatrix(degrees_in, degrees_out);
f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 400 400]);
imshow(full(A_scalefree));
print(f_scalefree, '../Figures/A_scalefree.png', '-dpng', '-r300')

%% Reconstructing A:
% Using A = spalloc(N, N, nonzeros); takes 0.142176 seconds.
% Using A = sparse(xidx, yidx, ones(nonzeros, 1, 'logical')); takes 0.049687 seconds.

clc;
tic;
nonzeros = sum(degrees_in);
xidx = uint16(zeros(nonzeros, 1));
yidx = uint16(zeros(nonzeros, 1));
choosefrom = uint16(2:N);
prob_leftout = degrees_out(1);
probs = degrees_out(choosefrom);

start = 1;
for i = 1:N
    num = degrees_in(i);
    
    % Sample with replacement when enough probabilities are nonzero.
    try
        % sampling without replacement
        chosen = datasample(choosefrom, num, 'Replace', false, 'Weights', probs.^2);
    catch
        % ok, there's not enough probabilities available
        % chosen = maxk(choosefrom, num);
        chosen = datasample(choosefrom, num, 'Replace', false, 'Weights', probs+1);
    end

    xidx(start:(start+(num-1))) = i;
    yidx(start:(start+(num-1))) = chosen;
    
    [~, chosenidx] = ismember(chosen, choosefrom);
    uniqueidx = unique(chosenidx);
    probs(uniqueidx) = probs(uniqueidx) - 1;    
    
    if i ~= N
        choosefrom(i) = i;
        
        tmp = probs(i);
        probs(i) = prob_leftout;
        prob_leftout = tmp;
    end
    start = start + num; 
end
A = sparse(xidx, yidx, ones(nonzeros, 1, 'logical'));
assert(full(sum(diag(A))) == 0)
whos A
toc

assert(sum(diag(A)) == 0)


imshow(full(A));

diffrows = degrees_in' - full(sum(A,2))';
diffcols = degrees_out' - full(sum(A,1));

N2 = N^2; thresh = 1.0e-9;
assert(sum(diffrows)/N2 < thresh)
assert(sum(diffcols)/N2 < thresh)

