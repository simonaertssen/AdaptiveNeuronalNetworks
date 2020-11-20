clear all; close al; clc;
% In this script we will add the functionality to check whether a given
% adjacency matrix consists of one or more connected components: are there
% parts of the network that are not connected to the rest?
% For an undirected graph A, the multiplicity m of the eigenvalue 0 of the 
% laplacian L gives us the number of connected components. This is only
% true for undirected graphs, but we can take the OR operator on the lower
% and upper triangular part of A to model that there actually exists a
% connection.

%% Setup
addpath('../Functions');
addpath('../Mean Field Reductions');

p.N = 10;
sfpars = make_scalefreeparameters(p, 3);
A = full(adjacencymatrix(sfpars.degrees_i, sfpars.degrees_o))

%% Check for connectivity:
% A first simple check is to see whether any of the row or columns sums are
% zero:
assert(all(sum(A,1) ~= 0))
assert(all(sum(A,2) ~= 0))

% Now we can take the OR operation:


