function P = randcauchy(seed, mu, gamma, n, m)
rng(seed);
    switch nargin
        case 4
            m = 1;
        case 3
            n = 1; m = 1;
        otherwise
            n = 1; m = 1; gamma = 1;
    end
    nsamples = 1000;
    P = mu + gamma*tan(pi*(linspace(0, 1, nsamples) - 0.5));
    P = repmat(P + sqrt(nsamples/N)*rand(nsamples,1), 1, N/nsamples);
end

