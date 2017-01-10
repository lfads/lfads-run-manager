function vec = makecol( vec )

% transpose if it's currently a row vector (unless its 0 x 1, keep as is)
if (size(vec,2) > size(vec, 1) && isvector(vec)) && ~(size(vec, 1) == 0 && size(vec, 2) == 1)
    vec = vec';
end

if size(vec, 1) == 1 && size(vec, 2) == 0
    vec = vec';
end

end

