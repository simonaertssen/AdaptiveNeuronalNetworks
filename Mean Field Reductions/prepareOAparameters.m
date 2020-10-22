function p = prepareOAparameters(p)
    p.k = unique(p.degrees_i);
    p.l = numel(p.k);
    pkperm = p.k(randperm(p.l));
    p.OA = zeros(p.l, p.l);
    for i = 1:p.l
        ks = p.k(i)*ones(p.l,1);
        p.OA(i, :) = p.P(p.k).*assortativity(p.k, pkperm, ks, ks, p.N, p.meandegree, 0)/p.meandegree;
    end
end