function v = reverseDims(v)
    % switch from col major to row major
    if ~isvector(v)
        v = permute(v, ndims(v):-1:1);
    end
end