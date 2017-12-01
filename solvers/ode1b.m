function [t, y] = ode1b(func, t, y0, options)
%% Default Arguments
if nargin == 3
    options = optimoptions('fsolve','Display','off');
end


%% Solve - Backward Euler
y = zeros(length(y0),length(t));
y(:,1) = y0;
for it = 2:length(t)
    [y(:,it),~,exitflag] = fsolve(@(y_new) y_new - y(:,it-1) -...
        (t(it)-t(it-1)) * func(t(it),y_new), y(:,it-1), options);
    if exitflag <= 0
        break;
    end
end

end