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


N = 32;
[~, Cu] = BVP_AD_upwind(L, U, Gamma, CE, N);
[x, Cc] = BVP_AD_central(L, U, Gamma, CE, N);
xh = zeros(N+1, 1); % fluxes are defined inbetween the nodes
phic_a = zeros(N+1, 1); phic_d = zeros(N+1, 1);
phiu_a = zeros(N+1, 1); phiu_d = zeros(N+1, 1);
for i = 1:N+1
    xh(i) = 0.5 * (x(i) + x(i+1));
    phic_a(i) = 0.5 * U * (Cc(i) + Cc(i+1));
    phic_d(i) = -Gamma * (Cc(i+1) - Cc(i)) / (x(i+1)-x(i));
    phiu_a(i) = U * Cu(i);
    phiu_d(i) = - Gamma * (Cu(i+1) - Cu(i)) / (x(i+1)-x(i));
end
% exact fluxes
Cex = @(x) (exp(Pe*x/L)-1)/(exp(Pe)-1);
phia_ex = @(x) U * Cex(x);
phid_ex = @(x) -Gamma * Pe/L *exp(Pe*x/L)/(exp(Pe)-1);
% plot discrete and exact fluxes together
figure
subplot(2,1,1)
plot(xh, phic_a, 'ob', xh, phic_d, 'sr', ...
xf, phia_ex(xf), '-k', xf, phid_ex(xf), '-k')
xlabel('x'); ylabel('\phi')
ylim([-2 2])
title('central')
subplot(2,1,2)
plot(xh, phiu_a, 'ob', xh, phic_d, 'sr', ...
xf, phia_ex(xf), '-k', xf, phid_ex(xf), '-k')
ylim([-2 2])
xlabel('x'); ylabel('\phi')
legend('\phi_c', '\phi_d', 'exact', 'location', 'NW')
title('upwind')

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

function [x, C] = BVP_AD_upwind(L, U, Gamma, CE, N)
    dx = L/N;
    x = (-dx/2:dx:L+dx/2)';
    A = zeros(N+2);
    b = zeros(N+2, 1);
    % interior points
    r1 = U/dx;
    r2 = Gamma / dx^2;
    for i = 2:N+1
        A(i, i-1) = r2 + r1;
        A(i, i) = -2 * r2 - r1;
        A(i, i+1) = r2;
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
