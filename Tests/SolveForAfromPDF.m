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
% Try to just always get the largest k values to put in a randomly chosen
% row.

N = 10;
pars.N = N;
% networkpars = make_scalefreeparameters(pars, 3);
networkpars = make_randomparameters(pars, 0.2);

degrees_in = networkpars.degrees_in;
degrees_out = networkpars.degrees_out;

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

choosefrom = cast(2:N,numtype);
prob_leftout = degrees_out(1);
probs = degrees_out;

din = degrees_in'
dout = degrees_out'

A = zeros(N,N);

rowpermutation = randperm(N);
for i = 1:N
    rowindex = i %= rowpermutation(i)
    numelements = degrees_in(rowindex)
    if numelements == 0
        continue
    end
    % Take out the diagonal element
%     allowedindices = [1:rowindex-1,rowindex+1:N];
%     addonesomewhere = allowedindices(randperm(N-1, 1))
    disp(probs')
    prob_leftout = probs(rowindex);
    probs(rowindex) = -1;
    [v, maxindex] = max(probs);
    probs(maxindex) = probs(maxindex) + 1;
    disp(probs')
    
    probs = flip(probs);
    disp(probs')
    [chosen, chosenidx] = maxk(probs, numelements)
    probs = flip(probs);
    [~, chosenidx] = ismember(chosen, probs)
    disp(probs')
    
    indices = idxidx(rowindex):idxidx(rowindex+1)-1;
    xidx(indices) = rowindex;
    yidx(indices) = chosenidx;
    
%     [~, chosenidx] = ismember(chosen, choosefrom);
%     uniqueidx = unique(chosenidx);
    probs(chosenidx) = probs(chosenidx) - 1;
    
    % Reset the probability vector:
    
    probs(rowindex) = prob_leftout;
    probs(maxindex) = probs(maxindex) - 1;

%     A(nonzeros(xidx(indices)), nonzeros(yidx(indices))) = 1
    
    probs
end

A = sparse(xidx, yidx, ones(numnonzeros, 1, 'logical'));
C = cat(1, full(A), degrees_out');
C = cat(2, C, [degrees_in; -100])

assert(sum(diag(A)) == 0);

diffrows = degrees_in' - full(sum(A,2))'
diffcols = degrees_out' - full(sum(A,1))


%%
% Create the matrix as sparse
if gpuDeviceCount > 0
    A = zeros(N, N, 'double');
    A(sub2ind([N, N], xidx, yidx)) = 1;
else
    A = sparse(xidx, yidx, ones(numnonzeros, 1, 'logical'));
end
assert(sum(diag(A)) == 0);

diffrows = degrees_in' - full(sum(A,2))';
diffcols = degrees_out' - full(sum(A,1));

N2 = N^2; thresh = 1.0e-6;
if sum(diffrows)/N2 > thresh
    sprintf(['Adjacency matrix might not be accurate: residue ', num2str(sum(diffrows)/N2)])
end
if sum(diffcols)/N2 > thresh
    sprintf(['Adjacency matrix might not be accurate: residue ', num2str(sum(diffcols)/N2)])
end

%% Test a fixed degree network:
netdegree = 300;
degrees_in = zeros(N,1);
degrees_in(randperm(N)) = netdegree;
degrees_out = degrees_in(randperm(N)); % because these two have to contain the same number of elements

assert(sum(degrees_in) == sum(degrees_out))

A_fixeddegree = adjacencymatrix(degrees_in, degrees_out);
f_fixeddegree = figure('Renderer', 'painters', 'Position', [0 800 400 400]);
hAxes = axes(f_fixeddegree); 
imshow(full(A_fixeddegree), 'Parent', hAxes);
title(hAxes, ['Fixed degree $$A_{ij}$$: $$N$$ = ', num2str(N), ', $$\langle k \rangle$$ = ', num2str(mean(degrees_in))],'interpreter','latex', 'FontSize', 15)
print(f_fixeddegree, '../Figures/A_fixeddegree.png', '-dpng', '-r300')

close(f_fixeddegree)

%% Test using the poisson distribution of random networks
netp = 0.60143;
meandegree = netp*(N - 1);
numlinks = netp*N*(N-1)/2;

degrees_in = poissrnd(meandegree,N,1);
degrees_out = degrees_in(randperm(N)); 

assert(sum(degrees_in) == sum(degrees_out))

A_random = adjacencymatrix(degrees_in, degrees_out);
f_random = figure('Renderer', 'painters', 'Position', [50 800 400 400]);
hAxes = axes(f_random); 
imshow(full(A_random), 'Parent', hAxes);
title(hAxes, ['Random $$A_{ij}$$: $$N$$ = ', num2str(N), ', $$\langle k \rangle$$ = ', num2str(round(mean(degrees_in)))],'interpreter','latex', 'FontSize', 15)
print(f_random, '../Figures/A_random.png', '-dpng', '-r300')

close(f_random)

%% Now using scale free networks:
pars.N = N;
degree = 3;
scalefreepars = make_scalefreeparameters(pars, degree);

degrees_out = scalefreepars.degrees_in(randperm(N)); 
assert(sum(scalefreepars.degrees_in) == sum(degrees_out))

A_scalefree = adjacencymatrix(scalefreepars.degrees_in, degrees_out);
f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 400 400]);
hAxes = axes(f_scalefree); 
imshow(full(A_scalefree), 'Parent', hAxes);
title(hAxes, ['Scale free $$A_{ij}$$: $$N$$ = ', num2str(N), ', $$ k \in $$ [', num2str(scalefreepars.kmin), ',', num2str(scalefreepars.kmax), '], $$\gamma$$ = ', num2str(degree)],'interpreter','latex', 'FontSize', 15)
print(f_scalefree, '../Figures/A_scalefree.png', '-dpng', '-r300')

close(f_scalefree)

%% Reconstructing A:
% Using A = spalloc(N, N, nonzeros); takes 0.142176 seconds.
% Using A = sparse(xidx, yidx, ones(nonzeros, 1, 'logical')); takes 0.049687 seconds.

clc;
tic;
numnonzeros = sum(degrees_in);
xidx = uint16(zeros(numnonzeros, 1));
yidx = uint16(zeros(numnonzeros, 1));
choosefrom = uint16(2:N);
prob_leftout = degrees_out(1);
probs = degrees_out(choosefrom);

start = 1;
for i = 1:N
    numelements = degrees_in(i);
    
    % Sample with replacement when enough probabilities are nonzero.
    try
        % sampling without replacement
        chosen = datasample(choosefrom, numelements, 'Replace', false, 'Weights', probs.^2);
    catch
        % ok, there's not enough probabilities available
        % chosen = maxk(choosefrom, num);
        chosen = datasample(choosefrom, numelements, 'Replace', false, 'Weights', probs+1);
    end

    xidx(start:(start+(numelements-1))) = i;
    yidx(start:(start+(numelements-1))) = chosen;
    
    [~, chosenidx] = ismember(chosen, choosefrom);
    uniqueidx = unique(chosenidx);
    probs(uniqueidx) = probs(uniqueidx) - 1;    
    
    if i ~= N
        choosefrom(i) = i;
        
        tmp = probs(i);
        probs(i) = prob_leftout;
        prob_leftout = tmp;
    end
    start = start + numelements; 
end
A = sparse(xidx, yidx, ones(numnonzeros, 1, 'logical'));
assert(full(sum(diag(A))) == 0)
whos A
toc

assert(sum(diag(A)) == 0)


imshow(full(A));

diffcols = degrees_in' - full(sum(A,2))';
diffrows = degrees_out' - full(sum(A,1));

N2 = N^2; thresh = 1.0e-9;
assert(sum(diffcols)/N2 < thresh)
assert(sum(diffrows)/N2 < thresh)
