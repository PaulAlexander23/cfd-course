
CE = 20;
C0 = 5;
K = 0.04;

tS = 6;
tL = 144;
t = 0:tS:tL;

CExact = CE + (C0 - CE)*exp(-K*t);
[~, CRK4] = rk4(@(t, y) K*(CE - y),t,C0);
[~, AB2] = ab2(@(t, y) K*(CE - y),t,C0);
[~, AM2] = am2(@(t, y) K*(CE - y),t,C0);

plot(t, CExact, t, CRK4, t, AB2, t, AM2);

figure
c = plot_convergence(@rk4, 12*2.^(0:-1:-7)');
c(2)

function c = plot_convergence(method, tS)
    
    CE = 20;
    C0 = 5;
    K = 0.04;
    
    tL = 24;
    CExact = CE + (C0 - CE)*exp(-K*tL);
    error = zeros(length(tS),1);
    
    for tSI = 1:length(tS)
        t = 0:tS(tSI):tL;
        
        [~, CApprox] = method(@(t, y) K*(CE - y),t,C0);
        
        error(tSI) = abs(CApprox(end) - CExact);
        
    end
    scatter(log10(tS), log10(error));
    hold on
    X = [ones(length(tS),1), log10(tS)];
    c = X \ log10(error);
    plot(log10(tS), X*c);
end