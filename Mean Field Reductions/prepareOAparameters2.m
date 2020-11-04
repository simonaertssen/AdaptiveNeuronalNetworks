function p = prepareOAparameters2(p)

    [p.k, ~, ic] = unique([p.degrees_i, p.degrees_o], 'rows');
    p.kcount = accumarray(ic, 1);
    p.Mk = numel(p.k)/2;

    p.OA = zeros(p.Mk, p.Mk);
    for i = 1:p.Mk
        p.OA(i, :) = p.P(p.k(:,1)).*p.P(p.k(:,2)).*assortativity(p.k(:,1), p.k_o(:,1), p.k(i,1), p.k_o(i,2), p.N, p.meandegree, 0)*p.N;
    end
    p.OA = p.K*p.OA/p.meandegree;
end
