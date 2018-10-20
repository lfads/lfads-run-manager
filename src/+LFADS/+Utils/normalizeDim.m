function tensor = normalizeDim(tensor, dim)

    ssq = sum(tensor.^2, dim);
    tensor = tensor ./ sqrt(ssq);

end