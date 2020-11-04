function scalefreepars = make_scalefreeparameters(pars, degree, kmin, kmax)
    scalefreepars = pars;
    % Make the necessery parameters for the scalefree networks
    if nargin < 1; error('Not enough input arguments'); end
    if nargin < 2; degree = 3; end
    if nargin < 3; kmin = round(pars.N*3/20); end
    if nargin < 4; kmax = round(pars.N*2/5); end
    assert(kmin < kmax);
   
    if degree < 2
        error('Scale free networks do not exist for degress less than 2');
    end
    
    scalefreepars.kmin = kmin;
    scalefreepars.kmax = kmax;
    scalefreepars.degree = degree;
    
    scalefreepars.P = @(x) scalefreepdf(x, pars.N, degree, kmin, kmax);
%     normalisation = scalefreepars.N/sum(scalefreepars.P(kmin:kmax));
%     scalefreepars.P = @(x) scalefreepars.P(x)*normalisation;
    
    % Improve the support of P for de unique degree vector k: take kmin:kmax
    idx = randperm(pars.N); kminmax = kmin:kmax; n = numel(kminmax);
    scalefreepars.degrees_i = zeros(pars.N,1);
    scalefreepars.degrees_i(idx(1:n)) = kmin:kmax;
    scalefreepars.degrees_i(idx(n+1:end)) = randsample(kmin:kmax, pars.N-n, true, scalefreepars.P(kmin:kmax))';
    
    if max(scalefreepars.degrees_i) > pars.N - 1
        disp(['Setting higher in-degrees to ', num2str(pars.N-1)]);
        scalefreepars.degrees_i(scalefreepars.degrees_i > pars.N - 1) = pars.N - 1;
    end
    scalefreepars.degrees_o = scalefreepars.degrees_i(randperm(pars.N));
 
    fsolveoptions = optimset('Display','off');
    scalefreepars.meandegree = fsolve(@(z) scalefreepars.P(z) - mean(scalefreepars.P(kmin:kmax)), kmin, fsolveoptions);
%     scalefreepars.meandegree = fsolve(@(z) sum(scalefreepars.P(kmin:z)) - sum(scalefreepars.P(z+1:kmax)), kmin, fsolveoptions);
%     scalefreepars.meandegree = mean(scalefreepars.degrees_i);
end

