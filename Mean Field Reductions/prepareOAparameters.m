function p = prepareOAparameters(p)
    [p.k, ~, ic] = unique(p.degrees_i);
    p.kcount = accumarray(ic, 1);
    p.Mk = numel(p.k);
    p.k_o = unique(p.degrees_o);
    p.OA = zeros(p.Mk, p.Mk);
    for i = 1:p.Mk
%         p.OA(i, :) = p.P(p.k).*assortativity(p.k, p.k_o, p.k(i), p.k_o(i), p.N, p.meandegree, 0);
        p.OA(i, :) = p.P2D(p.k(i), p.k_o(i)).*assortativity(p.k, p.k_o, p.k(i), p.k_o(i), p.N, p.meandegree, 0);
    end
    p.OA = p.K*p.OA/p.meandegree;
end
