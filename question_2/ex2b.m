%% i

omega0 = 2*pi;
zeta = 0.1;
g = @(t) 0;

x0 = 1;
v0 = 1;

tN = 101;
tL = 10;
t = linspace(0, tL, tN);

omega = omega0 * sqrt(1 - zeta^2);
xExact = exp(-zeta*omega0*t) .* (x0 * cos(omega*t) + (zeta*x0 + v0/omega0)*...
    omega0^2/omega^2 * sin(omega*t));

[~, xApprox] = rk4(@(t, y) [y(2); -2*zeta*omega0*y(2) - omega0^2*y(1) + g(t)],...
    t,[x0;v0]);

plot(t, xExact, t, xApprox(1,:))

%% ii
a = 3;
omegaf = 6;

g = @(t) a*sin(omegaf * t);

xHat = a/sqrt((omegaf^2 - omega0^2)^2 + 4 * zeta^2 * omega0^2 * omegaf^2);
[~, xApprox] = rk4(@(t, y) [y(2); -2*zeta*omega0*y(2) - omega0^2*y(1) + g(t)],...
    t,[x0;v0]);

figure
plot(t, xApprox(1,:), [t(1),t(end)], [xHat, -xHat; xHat, -xHat])

%% iii

kappa = 1;

[~, xApprox] = rk4(@(t, y) [y(2); -kappa*abs(y(2))*y(2) - omega0^2*y(1) + g(t)],...
    t,[x0;v0]);

figure
plot(t, xApprox(1,:))
