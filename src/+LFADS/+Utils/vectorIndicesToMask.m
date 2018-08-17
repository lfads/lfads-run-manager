function mask = vectorIndicesToMask(idx, len)
    if islogical(idx)
        mask = LFADS.Utils.makecol(idx);
    else
        if ~isscalar(len) && isvector(len), len = numel(len); end % assume its a mask of the same size
        mask = false(len, 1);
        mask(idx) = true;
    end
end