function p = prepareOAparameters(p)
    if isfield(p, 'P2D')
        vec = min(p.degrees_i):2:max(p.degrees_i);
        [p.k,p.k_o] = meshgrid(vec, vec);
        p.k = reshape(p.k,[],1); p.k_o = reshape(p.k_o,[],1);
    else
        [p.k, ~, ic] = unique(p.degrees_i);
        p.kcount = accumarray(ic, 1);
        p.k_o = unique(p.degrees_o);
    end
    p.Mk = numel(p.k);
    p.OA = zeros(p.Mk, p.Mk);
    
    if isfield(p, 'P2D')
    for i = 1:p.Mk
        p.OA(i, :) = p.P2D(p.k, p.k_o).*assortativity(p.k, p.k_o, p.k(i), p.k_o(i), p.N, p.meandegree, 0);
    end
    else
    for i = 1:p.Mk
        p.OA(i, :) = p.P(p.k).*assortativity(p.k, p.k_o, p.k(i), p.k_o(i), p.N, p.meandegree, 0);
    end
    end
    p.OA = p.K*p.OA/p.meandegree;
end
