function P = matlab_randcauchy(seed, x0, gamma, n, m)
rng(seed);
    switch nargin
        case 4
            m = 1;
        case 3
            n = 100; m = 1;
        otherwise
            n = 100; m = 1; gamma = 1;
    end
    pd = makedist('tLocationScale','mu',x0,'sigma',gamma,'nu',1);
    P = random(pd, n, m);
end

