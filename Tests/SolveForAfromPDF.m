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
% networkpars = make_scalefreeparameters(pars, 3);
networkpars = make_randomparameters(pars, 0.3);

degrees_in = networkpars.degrees_in;
degrees_out = networkpars.degrees_out;

% degrees_in = linspace(0,9,N)';
% degrees_out = linspace(0,9,N)';

%% The algorithm:
clc
N = numel(degrees_in);
if max(degrees_in) >= N
    error('Degree too large');
end
numnonzeros = sum(degrees_in);

% Test for laptop version or other:
if version('-release') == "2020a"
    numtype = 'uint16';
else
    numtype = 'double';
end

xidx = zeros(numnonzeros, 1, numtype);
yidx = zeros(numnonzeros, 1, numtype);
idxidx = cumsum([1; degrees_in]); % For indexing the idx and yidx vector
disp(idxidx')

% Adjust for zeros in first and last row
prob_leftout = degrees_out(1);
probs = degrees_out;

din = degrees_in'
dout = degrees_out'

A = zeros(N,N);

rowpermutation = randperm(N);
for i = 1:N
    rowindex = rowpermutation(i)
    numelements = degrees_in(rowindex)
    if numelements == 0
        continue
    end
    prob_leftout = probs(rowindex); % Take out diagonal element
    probs(rowindex) = -1;
    
    % Permutation makes the implementation quite robust: 
    % Don't just sample the first maximum elements 
    probsperm = randperm(N);
    [~, probsperminv] = sort(probsperm)
    disp('Shuffled')
    disp(probs(probsperm)')
    [chosen, chosenidx] = maxk(probs(probsperm), numelements)
    chosenidx = probsperm(chosenidx);
    
    indices = idxidx(rowindex):idxidx(rowindex+1)-1;
    xidx(indices) = rowindex;
    yidx(indices) = chosenidx;
    
    probs(chosenidx) = probs(chosenidx) - 1;
    
    % Reset the probability vector:
    probs(rowindex) = prob_leftout;

%     A(nonzeros(xidx(indices)), nonzeros(yidx(indices))) = 1
end

A = sparse(xidx, yidx, ones(numnonzeros, 1, 'logical'));
A(N,N) = 0;
assert(sum(diag(A)) == 0);

diffcols = degrees_out' - full(sum(A,1));
nonzeroidx = find(diffcols);
if numel(nonzeroidx) > 0
    
end


% C = cat(1, full(A), degrees_out');
% C = cat(2, C, [degrees_in; -100])
% 
% assert(sum(diag(A)) == 0);
% 
% diffrows = degrees_in' - full(sum(A,2))'
% diffcols = degrees_out' - full(sum(A,1))

%% Test the function:

