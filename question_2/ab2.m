function [t, y] = ab2(f, t, y0)
    
    y = zeros(length(y0),length(t));
    y(:,1) = y0;
    y(:,2) = y(:,1) + (t(2)-t(1))*f(t(1), y(:,1));
    
    for tI = 2:length(t)-1
        tS = t(tI+1)-t(tI);
        y(:,tI+1) = y(:,tI) +...
            tS/2 * (3*f(t(tI), y(tI)) - f(t(tI-1), y(tI-1)));
    end
    
end