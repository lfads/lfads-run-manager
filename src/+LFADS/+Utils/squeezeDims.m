function tsq = squeezeDims(t, dims)
    % like squeeze, except only collapses singleton dimensions in list dims
    siz = size(t);
    dims = dims(dims <= ndims(t));
    dims = dims(siz(dims) == 1);
    siz(dims) = []; % Remove singleton dimensions.
    siz = [siz ones(1,2-length(siz))]; % Make sure siz is at least 2-D
    tsq = reshape(t,siz);
end