function out = bsxfun_object(fn, a, b)
    % like bsxfun, except works with any type compatible with repmat

    % build size vectors with same length
    sza = size(a);
    szb = size(b);
    D = max(numel(sza), numel(szb));
    sza = cat(2, sza, ones(1, D-numel(sza)));
    szb = cat(2, szb, ones(1, D-numel(szb)));

    if any(sza ~= 1 & szb ~= 1 & sza ~= szb)
        error('Non-singleton dimension sizes must match');
    end

    repa = ones(1, D);
    repb = ones(1, D);

    mask = sza == 1 & szb ~= 1;
    repa(mask) = szb(mask);

    mask = sza ~= 1 & szb == 1;
    repb(mask) = sza(mask);

    a = repmat(a, repa);
    b = repmat(b, repb);

    out = arrayfun(fn, a, b);

end
