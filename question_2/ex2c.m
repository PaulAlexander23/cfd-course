
f0 = [0; 0; 0.4696];

etaL = 6;

[eta, f] = ode45(@(t, y) [y(2);y(3);-y(1)*y(3)], [0, etaL], f0);

plot(eta,f)

%% iii

Uinfty = 0.1;
nu = 1e-6;

x = linspace(0,100, 200);
z = linspace(0,0.4);
[X, Z] = meshgrid(x,z);

Eta = Z .* sqrt(Uinfty ./ (2* nu * X));
Eta(Eta > etaL) = etaL;
u = interp1(eta, Uinfty * f(:,2), Eta);

figure
contourf(X, Z, u);
xlabel('x [m]');
xlabel('z [m]');
title('Contours of velocity')
colorbar