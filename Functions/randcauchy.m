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
    P = reshape(mu + gamma*tan(pi*(linspace(1.e-6, 1-1.e-6, n) - 0.5)), n, m);
end

