function p = prepareOAparameters(p)
    p.k = unique(p.degrees_i);
    p.l = numel(p.k);
    pkout = unique(p.degrees_o);
    p.OA = zeros(p.l, p.l);
    for i = 1:p.l
        ks = p.k(i)*ones(p.l,1);
%         p.OA(i, :) = p.P(pkout).*assortativity(p.k, pkout, ks, pkout, p.N, p.meandegree, 0)/p.meandegree;
        p.OA(i, :) = p.P(pkout).*assortativity(p.k, pkout, ks, pkout, p.N, p.meandegree, 0)/p.meandegree;
    end
end