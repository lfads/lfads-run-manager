function map = hslmap(n, varargin)
% map = hslmap(n, ...)
% parameters: fracHueShift, saturation, luminance, fracHueSpan
% returns hue-spaced color map in rgb
    
    p = inputParser();
    p.addParameter('fracHueShift', -0.02, @isscalar); % between 0 and 1
    p.addParameter('saturation', 0.95, @isscalar);
    p.addParameter('luminance', 0.65, @isscalar);
    p.addParameter('fracHueSpan', 0.9, @isscalar);
    p.parse(varargin{:});
    
    hues = mod(p.Results.fracHueShift + circspace(0, p.Results.fracHueSpan, n)', 1);    
    hsl = [hues, p.Results.saturation * ones(n, 1), p.Results.luminance * ones(n,1)];
    map = hsl2rgb(hsl);
end

% adapted from colorspace.m utility by Pascal Getreuer 2005-2010
function out = hsl2rgb(map)
 % Convert HSL to sRGB
   L = map(:,3);
   Delta = map(:,2).*min(L,1-L);
   out = huetorgb(L-Delta,L+Delta,map(:,1)*360);
end

function Image = huetorgb(m0,m2,H)
    % Convert HSV or HSL hue to RGB
    N = size(H, 1);
    H = min(max(H(:),0),360)/60;
    m0 = m0(:);
    m2 = m2(:);
    F = H - round(H/2)*2;
    M = [m0, m0 + (m2-m0).*abs(F), m2];
    Num = length(m0);
    j = [2 1 0;1 2 0;0 2 1;0 1 2;1 0 2;2 0 1;2 1 0]*Num;
    k = floor(H) + 1;
    Image = reshape([M(j(k,1)+(1:Num).'),M(j(k,2)+(1:Num).'),M(j(k,3)+(1:Num).')],[N,3]);
end

function v = circspace(d1, d2, n)
% v = circspace(d1, d2, n)
% like linspace, except considers d1 == d2 in a circular axis

    if nargin == 2
        n = 100;
    else
        n = floor(double(n));
    end

    delta = (d2-d1)/n;
    v = linspace(d1, d2-delta, n);

end