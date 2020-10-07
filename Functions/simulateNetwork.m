function [tout, xout] = simulateNetwork(makeA, on_gpu)
    if nargin < 3
        on_gpu = @arrayinit
    end
    
    A = makeA(on_gpu)
    F = @(t, x, e, KdivN, a) thetaneurons_adjacency(t, x, e, KdivN, a, A_fixeddegree);
    [tout, xout] = DOPRI_threshold(F, tnow, tend, IC, h, pars, on_gpu);

end

