function [t, y] = ode2b(func, t, y0, options)
%% Default Arguments
if nargin == 3
    options = optimoptions('fsolve','Display','off');
end

%% Solve - Crank Nicolson
y = zeros(length(y0),length(t));
y(:,1) = y0;
for it = 2:length(t)
    F = func(t(it),y(:,it-1));
    [y(:,it),~,exitflag] = fsolve(@(y_new) y_new - y(:,it-1) -...
        0.5 * (t(it)-t(it-1)) * (func(t(it),y_new) + F),...
        y(:,it-1), options);
    if exitflag <= 0
        break;
    end
end

end