function lognormpars = make_lognormparameters(pars, mu, sigma, kmin)
    % Make the necessery parameters for the random networks
    lognormpars = pars;
    if nargin < 1; error('Not enough input arguments'); end
    if nargin < 2; mu = 10; end
    if nargin < 3; sigma = sqrt(mu); end
    if nargin < 4; kmin = round(pars.N/2); end
    
    mean = exp(mu + sigma*sigma/2);
    lognormpars.kmin = kmin;
    lognormpars.degrees_i = kmin + lognrnd(mu, sigma, [pars.N,1]);
    % Set higher degrees to N - 1:
    if max(lognormpars.degrees_i) > pars.N - 1
        disp(['Setting higher degrees to ', num2str(pars.N-1)]);
        lognormpars.degrees_i(lognormpars.degrees_i > pars.N - 1) = pars.N - 1;
    end
    lognormpars.degrees_o = lognormpars.degrees_i(randperm(pars.N));
    normalisation = pars.N/sum(lognpdf(0:pars.N, mu, sigma));
    lognormpars.P = @(x) lognpdf(x-kmin, mu, sigma)*normalisation;  
    lognormpars.meandegree = kmin + mean;
end

