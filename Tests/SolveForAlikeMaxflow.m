close all; clear all; clc
addpath('../Functions');

%% Making the adjacency matrix from the maximum flow
% From P(k) ~ k_i and P(k) ~ k_j to A_ij: a combinatorial problem...
% k_in = sum(A,2), k_out = sum(A,1)
% A -> k_in
% |
% v
% k_out

% Here we will make a sink and a source and compute A from the maximum flow
% in the graph: source -> k_in -> nodes -> A -> nodes -> k_out -> sink

N = 20;
netdegree = round(0.3*N);
degrees_in = zeros(N,1);
degrees_in(randperm(N)) = netdegree;
degrees_out = degrees_in(randperm(N));

%% Make a graph with sources and sinks:
idx = uint16((1:N)+2);
[xidx,yidx] = meshgrid(idx, idx);
xidx(1:N+1:end)=[]; yidx(1:N+1:end)=[];

s = [ones(1,N,'uint16'), idx, xidx];
t = [idx, 2*ones(1,N,'uint16'), yidx];
weights = [degrees_out; degrees_in; ones(N*(N-1),1)];
G = digraph(s,t,weights);
[mf,GF,cs,ct] = maxflow(G,1,2);

figure; hold on;
H = plot(G,'Layout','force');
H.EdgeLabel = {};
highlight(H,GF,'EdgeColor','r','LineWidth',2);
st = GF.Edges.EndNodes;
labeledge(H,st(:,1),st(:,2),GF.Edges.Weight);

figure;
H = plot(G,'Layout','layered','Sources',cs,'Sinks',ct, ...
    'EdgeLabel',G.Edges.Weight);
highlight(H,cs,'NodeColor','red')
highlight(H,ct,'NodeColor','green')
%% Make a digraph:
idx = uint16(1:N);
[xidx,yidx] = meshgrid(idx, idx);
S = sparse(xidx,yidx,ones(N*N,1,'logical')) - speye(N);
G = digraph(S);

[mf,GF] = maxflow(G,1,5);
plot(G,'EdgeLabel',G.Edges.Weight,'Layout','layered');
H.EdgeLabel = {};
highlight(H,GF,'EdgeColor','r','LineWidth',2);
st = GF.Edges.EndNodes;
labeledge(H,st(:,1),st(:,2),GF.Edges.Weight);

