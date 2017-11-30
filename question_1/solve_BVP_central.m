function [C, x] = solve_BVP_central(N,L,Gamma,U, CE)
    %SOLVE_BVP_CENTRAL Summary of this function goes here
    %   Detailed explanation goes here
    A = zeros(N+2);
    b = zeros(N+2, 1);
    
    dx = L/N;
    x = -dx/2:dx:L+dx/2;
    
    % interior points
    r = Gamma / dx^2;
    u = U / dx;
    for i = 2:N+1
        A(i, i-1) = u/2 + r;
        A(i, i)   = -2 * r;
        A(i, i+1) = -u/2 + r;
        b(i) = 0;
    end
    
    % BC1
    A(1,1) = 0.5;
    A(1,2) = 0.5;
    b(1) = 0;
    
    % BC2
    A(N+2,N+1) = 0.5;
    A(N+2,N+2) = 0.5;
    b(N+2) = CE; % BC2
    
    C = A \ b; % invert matrix
end

