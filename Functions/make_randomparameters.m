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
    randompars.degrees_i = poissrnd(randompars.meandegree, [pars.N,1]);
    % Set higher degrees to N - 1:
    if max(randompars.degrees_in) > pars.N - 1
        disp(['Setting higher degrees to ', num2str(pars.N-1)]);
        randompars.degrees_i(randompars.degrees_in > pars.N - 1) = pars.N - 1;
    end
    randompars.degrees_o = randompars.degrees_i(randperm(pars.N));
    randompars.P = @(x) poisspdf(x, randompars.meandegree)*pars.N;   
end

