function [xn, intervals] = bisection(f,xl,xr,tau,Nmax)
% Solves f(x) = 0 via bisection
% Inputs: 
%    f     function to find root of
%    xl,xr initial search interval bounds
%    tau   tolerance for convergence
%    Nmax  max number of iterations
% Outputs:
%    xn         sequence of iterates
%    intervals  sequence of bisection intervals (for visualization only)

% will return an error if initial interval condition is not met
assert(sign(f(xl))~=sign(f(xr)),'Interval ends do not have opposite signs.')  

xn = zeros(1,Nmax);        % pre allocate space for iterates
intervals = zeros(2,Nmax);
intervals(:,1) = [xl; xr];

for i = 1:Nmax
    xn(:,i) = 0.5*(xl+xr);
    if norm(f(xn(:,i))) < tau
        break
    elseif sign(f(xl))==sign(f(xn(:,i)))
        xl = xn(:,i);
    else
        xr = xn(:,i);
    end
    
    intervals(:,i+1) = [xl; xr];
end
xn = xn(:,1:i);
intervals = intervals(:,1:i);
end