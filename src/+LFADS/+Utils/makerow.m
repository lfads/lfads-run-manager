function vec = makerow( vec )
% convert vector to row vector

% leave size == [1 0] alone too
if(size(vec, 1) > size(vec,2) && isvector(vec)) && ~(size(vec, 2) == 0 && size(vec, 1) == 1)
    vec = vec';
end

if size(vec, 1) == 0 && size(vec, 2) == 1
    vec = vec';
end

end

