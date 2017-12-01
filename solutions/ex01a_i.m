clear all;
close all;
% parameters
L = 1;
U = 1;
Pe = 10;
CE = 1;
Gamma = 1/Pe;
% exact solution
Cex = @(x) (exp(Pe*x/L)-1)/(exp(Pe)-1);
xf = linspace(0, L, 1000);

Ns = [4, 8, 16, 32];
figure
for m = 1:length(Ns)
    [x, Cc] = BVP_AD_central(L, U, Gamma, CE, Ns(m));
    dx = x(2)-x(1);
    subplot(2,2,m)
    plot(x, Cc, '-^b', xf, Cex(xf), '-k');
    xlabel('x/L')
    ylabel('C/C_E')
    if (m==1), legend('central', 'exact', 'location', 'NW'), end
    title(['N = ', num2str(Ns(m)), '; Pe_C=', num2str(U*dx/Gamma)])
    xlim([0 L])
end

function [x, C] = BVP_AD_central(L, U, Gamma, CE, N)
    dx = L/N;
    x = (-dx/2:dx:L+dx/2)';
    A = zeros(N+2);
    b = zeros(N+2, 1);
    % interior points
    r1 = 0.5*U/dx;
    r2 = Gamma / dx^2;
    for i = 2:N+1
        A(i, i-1) = r2 + r1;
        A(i, i) = -2 * r2;
        A(i, i+1) = r2 - r1;
        b(i) = 0;
    end
    % BC1
    A(1,1) = 0.5;
    A(1,2) = 0.5;
    b(1) = 0;
    % BC2
    A(N+2,N+1) = 0.5;
    A(N+2,N+2) = 0.5;
    b(N+2) = CE;
    C = A \ b; % invert matrix
end
