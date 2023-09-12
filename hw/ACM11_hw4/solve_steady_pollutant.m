function [u_mesh, ux_mesh] = solve_steady_pollutant(x_mesh, a, D, p)
    dx = x_mesh(2)- x_mesh(1);
    
    rhs = p(x_mesh');
    A = getMatrix(x_mesh,a,D);
    
    u_mesh = -A\rhs;
    ux_mesh = ([u_mesh(2); u_mesh(2:end)] - [u_mesh(1); u_mesh(1:end-1)])/dx;
end

function A = getMatrix(x_mesh,a,D)
    dx = x_mesh(2)- x_mesh(1);
    N = length(x_mesh);
    A1 = 1/2*diag(ones(N-1,1),1)- 1/2*diag(ones(N-1,1),-1);
    A2 = -2*diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
    A = - a/dx * A1 + D/(dx^2) * A2;
end