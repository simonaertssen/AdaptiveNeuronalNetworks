function p = prepareOAparameters_sf(p)
    %     [p.k, ~, ic] = unique(p.degrees_i);
    p.k = [p.kmin:p.kmax]';
%     p.kcount = accumarray(ic, 1);
    p.Mk = numel(p.k)

    p.k_o = p.k;
    p.OA = zeros(p.Mk, p.Mk);
    for i = 1:p.Mk
        p.OA(i, :) = p.P(p.k).*assortativity(p.k, p.k_o, p.k(i), p.k_o(i), p.N, p.meandegree, 0);
    end
    p.OA = p.K*p.OA/p.meandegree;
end
