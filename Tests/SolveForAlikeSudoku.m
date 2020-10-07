close all; clear all; clc

%% Making the adjacency matrix from a degree distribution
% From P(k) ~ k_i and P(k) ~ k_j to A_ij: a combinatorial problem...
% k_in = sum(A,2), k_out = sum(A,1)
% A -> k_in
% |
% v
% k_out
% The challenge is to get k_in(i) for each row, and then to sample from
% k_out so that we get the right number of degrees.

%% Test using the poisson distribution of random networks
N = 10;
netp = 0.3;
meandegree = netp*(N - 1);
numlinks = netp*N*(N-1)/2;

degrees_in = poissrnd(meandegree,N,1)';
degrees_out = degrees_in(randperm(N)); % because these two have to contain the same number of elements

P = @(x) poisspdf(x, meandegree)*N;   

%% Reconstructing A
% Using A = spalloc(N, N, nonzeros); takes 0.142176 seconds.
% Using A = sparse(xidx, yidx, ones(nonzeros, 1, 'logical')); takes 0.049687 seconds.
clc;
tic;
nonzeros = sum(degrees_in);
xidx = zeros(nonzeros, 1, 'single');
yidx = zeros(nonzeros, 1, 'single');
choosefrom = 2:N;
prob_leftout = degrees_out(1);
probs = degrees_out(choosefrom);

deleteme = zeros(N, N);
start = 1;
for i = 1:N
    num = degrees_in(i);
    xidx(start:(start+(num-1))) = i;
    
    % Sample with replacement when enough probabilities are nonzero.
    % This can only be done with a loop, one by one.
%     chosen = zeros(num,1);
%     for n = 1:num
%         choice = randsample(choosefrom,1,true,probs);
%         probs(choosefrom==choice) = probs(choosefrom==choice) - 1;
%         chosen(n) = choice;
%     end
    chosen = datasample(choosefrom, num, 'Replace', num > nnz(probs), 'Weights', probs)
%         choosefrom

    yidx(start:(start+(num-1))) = chosen;

    [~, chosenidx] = ismember(chosen, choosefrom)
    uniqueidx = unique(chosenidx)
    freq = unique(sum(uniqueidx==uniqueidx'), 'stable')
%     [freq, uniqueidx] = hist(chosenidx,unique(chosenidx))
%     frequency = sum(chosenidx==chosenidx')
%     [freq,uniqueidx] = groupcounts(chosenidx) 
    probs
    probs(uniqueidx) = probs(uniqueidx) - freq;
    probs
    deleteme(i*ones(num,1), choosefrom(randperm(N-1,num))) = 1;
    
    
    if i ~= N
        choosefrom(i) = i;
        
        tmp = probs(i);
        probs(i) = prob_leftout;
        prob_leftout = tmp;
    end
%     probs(probs < 0) = 0
    start = start + num; 
end
A = sparse(xidx, yidx, ones(nonzeros, 1, 'logical'));
whos A
toc


disp(degrees_in - full(sum(A,2))')
disp(degrees_out - full(sum(A,1)))

% immse(degrees_in, full(sum(A,2))')
% immse(degrees_out, full(sum(A,1))')

%Or go through and switch columns