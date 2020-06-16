function gr = symequidgrid(n, b)
% Create symmetrical equidistant grid. Supports one or more variables
% Parameters:
%   n   - #(s) of points
%   b   - bound(s) of variable(s) (just the positive one, assumed symmetrical around origin)

nn = length(n); nb = length(b);
if nb > 1 && nn == 1,
    n = n + zeros(nb, 1);
elseif nb ~= nn,
    error('n and b must have the same length, or n must be scalar');
end;

if nb == 1,
    gr = -b : 2*b/(n-1) : b;
else
    gr = cell(nb, 1);
    for i = 1:nb,
        if i == 1
        gr{i} = -b(i): b(i)/(n(i)-1) : b(i);
        end
        if i == 2
            gr{i} = 0;
        end
    end;
end;

end
