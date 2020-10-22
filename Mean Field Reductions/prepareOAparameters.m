function p = prepareOAparameters(p)
    p.k = unique(p.degrees_i, 'stable');
    p.l = numel(p.k);
    p.k_o = unique(p.degrees_o, 'stable');
    p.OA = zeros(p.l, p.l);
    for i = 1:p.l
        p.OA(i, :) = p.P(p.k).*assortativity(p.k, 0, 0, p.k_o(i), p.N, p.meandegree, 0)/p.meandegree;
        
%good:         p.OA(i, :) = p.P(p.k_o).*assortativity(p.k, p.k_o, ks, p.k_o, p.N, p.meandegree, 0)/p.meandegree;
%better:         p.OA(i, :) = p.P(p.k(i)).*assortativity(p.k(i), 0, 0, p.k_o, p.N, p.meandegree, 0)/p.meandegree;
% best:         p.OA(i, :) = p.P(p.k).*assortativity(p.k, 0, 0, p.k_o(i), p.N, p.meandegree, 0)/p.meandegree;

    end
end

% p.OA(i, :) = p.P(p.k).*assortativity(p.k_o, 0, 0, p.k_o(i), p.N, p.meandegree, 0)/p.meandegree;
