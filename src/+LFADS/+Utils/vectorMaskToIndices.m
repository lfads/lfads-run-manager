function idx = vectorMaskToIndices(mask)
    if islogical(mask)
        idx = LFADS.Utils.makecol(find(mask)); %#ok<FNDSB>
    else
        idx = LFADS.Utils.makecol(mask);
    end
end