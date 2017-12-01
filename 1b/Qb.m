

a = 1;
Cb = 1;
L = 20;

N = 2^4;
dx = L/N;

eta = (-dx/2:dx:L+dx/2)';

D1 = spdiags(ones(N+2,1)*[-0.5,0,0.5],[-1,0,1],N+2,N+2)/dx;
D2 = spdiags(ones(N+2,1)*[1,-2,1],[-1,0,1],N+2,N+2)/dx/dx;

%A = D1*((1+eta.^3).*D1);
A = (1+eta.^3).*D2 + 3*eta.^2.*D1;
b = zeros(N+2,1);

A(1,1) = -1/dx - a*0.5;
A(1,2) = 1/dx - a*0.5;
A(1,3) = 0;
b(1) = 0;

A(end,end-2) = 0;
A(end,end-1) = 0.5;
A(end,end) = 0.5;
b(end) = Cb;

C = A \ b;

eta2 = linspace(0,20);
f = @(eta) sqrt(3)/2/pi*(log((eta+1)./(sqrt(eta.^2 - eta + 1))) - sqrt(3)*(pi/2 - atan((2*eta-1)/sqrt(3))));
%Cex = Cb + (Cb - Cw)*f(eta2);

plot(eta,C,'-ob')
%plot(eta2,Cex,'k')
axis([0,20,-0.5,1])