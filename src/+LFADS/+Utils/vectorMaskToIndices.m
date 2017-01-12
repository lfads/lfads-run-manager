function idx = vectorMaskToIndices(mask)
    if islogical(mask)
        idx = makecol(find(mask));
    else
        idx = makecol(mask);
    end
end