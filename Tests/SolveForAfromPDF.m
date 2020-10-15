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

% Adjust for zeros in first and last row
prob_leftout = degrees_out(1);
probs = degrees_out;

din = degrees_in'
dout = degrees_out'

A = zeros(N,N);

rowpermutation = randperm(N);
for i = 1:N
    rowindex = i;
    numelements = degrees_in(rowindex);
    if numelements == 0
        continue
    end
    prob_leftout = probs(rowindex); % Take out diagonal element
    probs(rowindex) = -1;
    
    % Permutation makes the implementation quite robust: 
    % Don't just sample the first maximum elements 
    probsperm = randperm(N);
    [~, probsperminv] = sort(probsperm);
    [chosen, chosenidx] = maxk(probs(probsperm), numelements);
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

C = cat(1, full(A), degrees_out');
C = cat(2, C, [degrees_in; -100])

assert(sum(diag(A)) == 0);

diffcols = degrees_out' - full(sum(A,1))
nonzeroidx = find(diffcols)
numel(nonzeroidx)
if numel(nonzeroidx) > 0
    if diffcols(nonzeroidx(1)) == -diffcols(nonzeroidx(2))
        switchidx1 = find(A(:,nonzeroidx(1)) == 1)
        for sw = 1:numel(switchidx1)
            sw
            A(sw,nonzeroidx(1)) 
            A(sw,nonzeroidx(2)) 
            if A(sw,nonzeroidx(2)) == 0
                A(sw,nonzeroidx(2)) = 1;
                A(sw,nonzeroidx(1)) = 0;
                break
            end
        end
    end
    C = cat(1, full(A), degrees_out');
    C = cat(2, C, [degrees_in; -100])
end

diffcols = degrees_out' - full(sum(A,1))

%% Test the function:
pars.N = 1000;
netp = 0.60143;
meandegree = netp*(N - 1);
networkpars = make_randomparameters(pars, netp);

assert(sum(networkpars.degrees_in) == sum(networkpars.degrees_out))

tic 
A_random = adjacencymatrix(networkpars.degrees_in, networkpars.degrees_out);
toc

tic 
A_random = adjacencymatrix_from_sampling(networkpars.degrees_in, networkpars.degrees_out);
toc

% Lesson: faster if we just get it after one try, but it seems to be more
% accurate.

%% Test a fixed degree network:
pars.N = 500;
netdegree = 50;
fdpars = make_fixeddegreeparameters(pars, netdegree); 

A_fixeddegree = adjacencymatrix(fdpars.degrees_in, fdpars.degrees_out);
f_fixeddegree = figure('Renderer', 'painters', 'Position', [0 800 400 400]);
hAxes = axes(f_fixeddegree); 
imshow(full(A_fixeddegree), 'Parent', hAxes);
title(hAxes, ['Fixed degree $$A_{ij}$$: $$N$$ = ', num2str(pars.N), ', $$\langle k \rangle$$ = ', num2str(fdpars.meandegree)],'interpreter','latex', 'FontSize', 10)
exportpdf(f_fixeddegree, '../Figures/A_fixeddegree1.pdf', true);

close(f_fixeddegree)

netdegree = 300;
fdpars = make_fixeddegreeparameters(pars, netdegree); 

A_fixeddegree = adjacencymatrix(fdpars.degrees_in, fdpars.degrees_out);
f_fixeddegree = figure('Renderer', 'painters', 'Position', [0 800 400 400]);
hAxes = axes(f_fixeddegree); 
imshow(full(A_fixeddegree), 'Parent', hAxes);
title(hAxes, ['Fixed degree $$A_{ij}$$: $$N$$ = ', num2str(pars.N), ', $$\langle k \rangle$$ = ', num2str(fdpars.meandegree)],'interpreter','latex', 'FontSize', 10)
exportpdf(f_fixeddegree, '../Figures/A_fixeddegree2.pdf', true);

close(f_fixeddegree)


%% Test using the poisson distribution of random networks
netp = 0.10043;
rdpars = make_randomparameters(pars, netp);

A_random = adjacencymatrix(rdpars.degrees_in, rdpars.degrees_out);
f_random = figure('Renderer', 'painters', 'Position', [50 800 400 400]);
hAxes = axes(f_random); 
imshow(full(A_random), 'Parent', hAxes);
title(hAxes, ['Random $$A_{ij}$$: $$N$$ = ', num2str(pars.N), ', $$\langle k \rangle$$ = ', num2str(round(rdpars.meandegree))],'interpreter','latex', 'FontSize', 10)
exportpdf(f_random, '../Figures/A_random1.pdf', true);

close(f_random)

netp = netdegree/(pars.N-1);
rdpars = make_randomparameters(pars, netp);

A_random = adjacencymatrix(rdpars.degrees_in, rdpars.degrees_out);
f_random = figure('Renderer', 'painters', 'Position', [50 800 400 400]);
hAxes = axes(f_random); 
imshow(full(A_random), 'Parent', hAxes);
title(hAxes, ['Random $$A_{ij}$$: $$N$$ = ', num2str(pars.N), ', $$\langle k \rangle$$ = ', num2str(round(rdpars.meandegree))],'interpreter','latex', 'FontSize', 10)
exportpdf(f_random, '../Figures/A_random2.pdf', true);
close(f_random)

%% Now using scale free networks:
degree = 2.1;
sfpars = make_scalefreeparameters(pars, degree, 50, 400);

A_scalefree = adjacencymatrix(sfpars.degrees_in, sfpars.degrees_out);
f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 400 400]);
hAxes = axes(f_scalefree); 
imshow(full(A_scalefree), 'Parent', hAxes);
title(hAxes, ['Scale free $$A_{ij}$$: $$N$$ = ', num2str(pars.N), ', $$ k \in $$ [', num2str(sfpars.kmin), ',', num2str(sfpars.kmax), '], $$\gamma$$ = ', num2str(degree)],'interpreter','latex', 'FontSize', 10)
exportpdf(f_scalefree, '../Figures/A_scalefree1.pdf', true);

close(f_scalefree)

degree = 10;
sfpars = make_scalefreeparameters(pars, degree, 50, pars.N);

A_scalefree = adjacencymatrix(sfpars.degrees_in, sfpars.degrees_out);
f_scalefree = figure('Renderer', 'painters', 'Position', [50 800 400 400]);
hAxes = axes(f_scalefree); 
imshow(full(A_scalefree), 'Parent', hAxes);
title(hAxes, ['Scale free $$A_{ij}$$: $$N$$ = ', num2str(pars.N), ', $$ k \in $$ [', num2str(sfpars.kmin), ',', num2str(sfpars.kmax), '], $$\gamma$$ = ', num2str(degree)],'interpreter','latex', 'FontSize', 10)
exportpdf(f_scalefree, '../Figures/A_scalefree2.pdf', true);
close(f_scalefree)
