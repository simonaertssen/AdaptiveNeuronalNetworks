function [struct] = prepareOApars(struct)
    struct.assortativity = @(kaccent, k, c) max(0, min(1, (kaccent.*k/(struct.meandegree*struct.N) + c*(1))));

    struct.k = unique(struct.degrees);
    struct.l = length(struct.k);
    struct.OA = zeros(struct.l, struct.l);
    for i = 1:struct.l
        ks = struct.k(i)*ones(struct.l,1);
        struct.OA(i, :) = struct.P(struct.k).*struct.assortativity(struct.k,ks,struct.c)/struct.meandegree;
    end
end

