function p = prepareOAparameters(p)
    p.k = unique(p.degrees_i, 'stable');
    p.l = numel(p.k);
    p.k_o = unique(p.degrees_o, 'stable');
    p.OA = zeros(p.l, p.l);
    for i = 1:p.l
        p.OA(i, :) = p.P(p.k).*assortativity(p.k, p.k_o, p.k(i), p.k_o(i), p.N, p.meandegree, 0);
    end
    p.OA = p.OA/p.meandegree
end
