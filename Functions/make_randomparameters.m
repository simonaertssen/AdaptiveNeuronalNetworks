function randompars = make_randomparameters(pars, netp)
    % Make the necessery parameters for the random networks
    randompars = pars;
    if netp > 1
        error('The network degree is too high');
    elseif netp < 0
        error('The network degree is negative');
    end    
    
    randompars.netp = netp;
    randompars.meandegree = netp*(pars.N - 1);
    
    stddev = sqrt(randompars.meandegree);
    base = round(randompars.meandegree + -1.966*stddev:1.966*stddev);
    
    idx = randperm(pars.N); n = numel(base);
    randompars.degrees_i = zeros(pars.N,1);
    randompars.degrees_i(idx(1:n)) = base;
    randompars.degrees_i(idx(n+1:end)) = poissrnd(randompars.meandegree, [pars.N-n,1]);
    
    % Set higher degrees to N - 1:
    if max(randompars.degrees_i) > pars.N - 1
        disp(['Setting higher degrees to ', num2str(pars.N-1)]);
        randompars.degrees_i(randompars.degrees_i > pars.N - 1) = pars.N - 1;
    end
    randompars.degrees_o = randompars.degrees_i(randperm(pars.N));
    randompars.P = @(x) poisspdf(x, randompars.meandegree)*pars.N;   
end

