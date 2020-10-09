function p = prepareOAparameters(p)
    ASSortativity = @(kaccent, k, c) max(0, min(1, (kaccent.*k/(p.meandegree*p.N) + c*(1))));
    p.k = unique(p.degrees_in);
    p.l = numel(p.k);
    p.OA = zeros(p.l, p.l);
    for i = 1:p.l
        ks = p.k(i)*ones(p.l,1);
        p.OA(i, :) = p.P(p.k).*ASSortativity(p.k,ks,0)/p.meandegree;
        %p.OA(i, :) = p.P(p.k).*assortativity(p.k, p.k(randperm(p.l)), ks, ks(randperm(p.l)), p.N, p.meandegree, 0)/p.meandegree;
    end
end

% p.assortativity = @(kaccent, k, c) max(0, min(1, (kaccent.*k/(p.meandegree*p.N) + c*(1))));
% 
%     p.k = unique(p.degrees);
%     p.l = length(p.k);
%     p.OA = zeros(p.l, p.l);
%     for i = 1:p.l
%         ks = p.k(i)*ones(p.l,1);
%         p.OA(i, :) = p.P(p.k).*p.assortativity(p.k,ks,p.c)/p.meandegree;
%     end
% end

