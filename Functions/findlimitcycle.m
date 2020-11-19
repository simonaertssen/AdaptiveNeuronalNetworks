function ts = findlimitcycle(vector)
    vector = detrend(vector, 1);
    ts = [nan, nan];
    t_fmean = find(diff(sign(vector))); 
    dvector = diff(vector);

    tstart = nan;
    for i = 3:4
        if sign(dvector(t_fmean(end))) == sign(dvector(t_fmean(end-i)))
            tstart = t_fmean(end-i);
        end  
    end
    ts = [tstart, t_fmean(end)];
end

