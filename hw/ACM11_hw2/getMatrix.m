function A = getMatrix(N, L)
    dx = L/N;
    A = diag(-2/dx^2 * [ones(1,N-1)]) + diag(1/dx^2 * [ones(1,N-2)], 1) + diag(1/dx^2 * [ones(1,N-2)], -1);
end