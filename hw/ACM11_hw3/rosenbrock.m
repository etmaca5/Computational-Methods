function [f,g,H] = rosenbrock(x,y)
% "nargout" is the number of output arguments requested
% if we call this function and ask for 1 output, it will just return f
% if we ask for 2 outputs, we will get f and g
% if we ask for 3 outputs, we get f, g, and H

f = (1-x).^2 + (y-x.^2).^2;

if nargout >= 2 
    g = [-2 * (1 - x) - 4 * x * (y - x.^2); 2 * (y - x.^2)]; 
end

if nargout == 3
    H = [2 - (4 * y) + 12 * (x.^2), -4 * x; -4 * x, 2];
end
