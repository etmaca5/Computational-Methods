function [opt_path, conv_path] = Rosenbrock_newton_fixed(x0, y0, alpha, epsilon, Nmax)
    xn = x0;
    yn = y0;
    opt_path = zeros(2, Nmax+1);
    conv_path = zeros(1, Nmax+1);
    for n = 1 : Nmax
        % call rosenbrock.m at (xn,yn) to get the current fn, gn, and Hn
        [fn, gn, Hn] = rosenbrock(xn, yn);
        
        % store fn in nth entry of conv_path, and store xn,yn in nth column
        % of opt path
        conv_path(n) = fn;
        opt_path(:,n) = [xn, yn];
        
        % check to see if norm of gradient is below epsilon tolerance - if
        % so, break the for loop.
        if norm(gn) < epsilon
            break;
        end
        
        % check if Newton update is a descent direction. If not, change the
        % sign of the Newton update
        update = Hn \ gn;
        if (gn).' * update < 0
            update = -1 * update;
        end
        
        % update xn and yn with Newton update
        xn = xn - (alpha * update(1));
        yn = yn - (alpha * update(2));
    end
    opt_path(:, n+1) = [xn; yn];
    conv_path(:, n+1) = rosenbrock(xn,yn);
    opt_path = opt_path(:, 1:n+1);
    conv_path = conv_path(:, 1:n+1);
end