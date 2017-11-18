clear all; close all;

Gamma = 0.005;
U = 0.05;
L = 1;
Pe = U*L/Gamma;

N = 2.^(4:10);
epsilon = zeros(7,1);
for j = 1:7
    dx = L/N(j);
    CE = 1;
    
    x = -dx/2:dx:L+dx/2;
    A = zeros(N(j)+2);
    b = zeros(N(j)+2, 1);
    
    % interior points
    r = Gamma / dx^2;
    u = U / dx;
    for i = 2:N(j)+1
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
    A(N(j)+2,N(j)+1) = 0.5;
    A(N(j)+2,N(j)+2) = 0.5;
    b(N(j)+2) = CE; % BC2
    
    C = A \ b; % invert matrix
    
    Cex = CE*(exp(Pe*x'/L)-1)/(exp(Pe)-1);
    epsilon(j) = 1/N(j) * sum(abs(Cex-C));
    
end

loglog(1./N, epsilon, '^b');
xlabel('\Delta x / L')
ylabel('\epsilon')