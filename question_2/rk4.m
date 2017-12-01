function [t, y] = rk4(f, t, y0)
    
    y = zeros(length(y0),length(t));
    y(:,1) = y0;
    
    for tI = 1:length(t)-1
        tS = t(tI+1)-t(tI);
        k1 = tS * f(t(tI), y(:, tI));
        k2 = tS * f(t(tI)+tS/2, y(:, tI) + k1/2);
        k3 = tS * f(t(tI)+tS/2, y(:, tI) + k2/2);
        k4 = tS * f(t(tI)+tS, y(:, tI) + k3);
        y(:, tI+1) = y(:, tI) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
    
end
