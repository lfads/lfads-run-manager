function vec = lindelta(from, delta, N)
% LINDELTA generates vector with specified size and spacing
%   vec = LINDELTA(from, delta, N) generates a column vector of N
%   linearly spaced points from from to delta.
%
%   See also LINSPACE, COLON
%

vec = from:delta:(from+(N-1)*delta);

end