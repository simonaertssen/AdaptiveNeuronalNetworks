close all; clear all; clc;
addpath('../Functions');

%% Investigate the properties of different network graphs

pars.N = 100;
if gpuDeviceCount > 0
    d = gpuDevice(gpuDeviceCount);
    disp(d)
end

%% Graph properties of a fixed degree network / diracnet: CORRECT
nsamples = 10;
netdegrees = linspace(log(pars.N)^3, pars.N, nsamples); % Use log to make sure k >>> log(N)

for i = 1:nsamples
    fixeddegreepars = make_fixeddegreeparameters(pars, round(netdegrees(i)));
    A_fixeddegree = initarr adjacencymatrix(fixeddegreepars.degrees_in);

end


