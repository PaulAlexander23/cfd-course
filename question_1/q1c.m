
CE = 1;
Gamma = 0.005;
U = 0.05;
L = 1;
Pe = U*L/Gamma;

for N = [4,8,16,32]
    
    [Cc,x] = solve_BVP_central(N,L,Gamma,U,CE);
    
    Cu = solve_BVP_upwind(N,L,Gamma,U,CE);
    
    xf = linspace(0, L, 1000);
    
    Cex = CE*(exp(Pe*xf/L)-1)/(exp(Pe)-1);
    
    figure()
    plot(x, Cc, '-^b', x, Cu, '-vr', xf, Cex, '-k')
    axis([0,1,-0.5*(N==4),1]);
end


N = 2.^(4:10);
lN = length(N);
epsilonC = zeros(1,lN);
epsilonU = zeros(1,lN);
for i = 1:lN
    
    [Cc,x] = solve_BVP_central(N(i),L,Gamma,U,CE);
    
    [Cu,x2] = solve_BVP_upwind(N(i),L,Gamma,U,CE);
    
    Cex = CE*(exp(Pe*x'/L)-1)/(exp(Pe)-1);
    
    epsilonC(i) = 1/N(i) * sum(abs(Cex-Cc));
    epsilonU(i) = 1/N(i) * sum(abs(Cex-Cu));
end

loglog(1./N, epsilonC, '-^b');
hold on;
loglog(1./N, epsilonU, '-vr');
xlabel('\Delta x / L')
ylabel('\epsilon')