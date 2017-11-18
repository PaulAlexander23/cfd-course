clear all; close all;

Gamma = 0.005;
L = 1;
N = 20;
dx = L/N;
CE = 1;

x = -dx/2:dx:L+dx/2;
A = zeros(N+2);
b = zeros(N+2, 1);

% interior points
r = Gamma / dx^2;
for i = 2:N+1
  A(i, i-1) = r;
  A(i, i)   = -2 * r;
  A(i, i+1) = r;
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

xf = linspace(0, L, 1000); Cex = xf/L*CE;
plot(x, C, 'ob', xf, Cex, '-k');
xlabel('x/L')
ylabel('C/C_E')
legend('numerical', 'exact', 'location', 'NW')