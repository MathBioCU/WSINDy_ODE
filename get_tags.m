function tags = get_tags(polys,trigs,n)
    tags = [];
    for p = 1:length(polys)
        monom_powers = partitionNk(polys(p),n);
        tags = [tags;monom_powers];
    end
    for k=1:length(trigs)
        trig_inds = [-trigs(k)*1i*eye(n);trigs(k)*1i*eye(n)];
        tags = [tags; trig_inds];
    end
end
