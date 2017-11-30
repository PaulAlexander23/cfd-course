function [t, y] = am2(f, t, y0)
    
    y = zeros(length(y0),length(t));
    y(:,1) = y0;
    
    options = optimoptions('fsolve', 'display', 'none');
    
    for tI = 1:length(t)-1
        tS = t(tI+1)-t(tI);
        y(:,tI+1) = fsolve(@(x) x - y(:,tI) -...
            tS/2 * (f(t(tI+1), x) + f(t(tI), y(tI))), y(:,tI), options);
    end
    
end