clear all;
close all;
% parameters
L = 1;
U = 1;
Pe = 10;
CE = 1;
Gamma = 1/Pe;
% exact solution
Cex = @(x) (exp(Pe*x/L)-1)/(exp(Pe)-1);
xf = linspace(0, L, 1000);

Ns = 10000;
% non-sparse version
tStart = tic;
[x, Cc] = BVP_AD_central(L, U, Gamma, CE, Ns);
tElapsed = toc(tStart);
fprintf('Time (non-sparse): %.3f\n',tElapsed)
figure
plot(x, Cc, '-^b', xf, Cex(xf), '-k');

% Sparse version 1
tStart = tic;
[x, Cc] = BVP_AD_central_sparse(L, U, Gamma, CE, Ns);
tElapsed = toc(tStart);
fprintf('Time (sparse 1): %.3f\n',tElapsed)
figure
plot(x, Cc, '-^b', xf, Cex(xf), '-k');

% Sparse version 2
tStart = tic;
[x, Cc] = BVP_AD_central_sparse2(L, U, Gamma, CE, Ns);
tElapsed = toc(tStart);
fprintf('Time (sparse 2): %.3f\n',tElapsed)
figure
plot(x, Cc, '-^b', xf, Cex(xf), '-k');

function [x, C] = BVP_AD_central(L, U, Gamma, CE, N)
    dx = L/N;
    x = (-dx/2:dx:L+dx/2)';
    A = zeros(N+2);
    b = zeros(N+2, 1);
    % interior points
    r1 = 0.5*U/dx;
    r2 = Gamma / dx^2;
    for i = 2:N+1
        A(i, i-1) = r2 + r1;
        A(i, i) = -2 * r2;
        A(i, i+1) = r2 - r1;
        b(i) = 0;
    end
    % BC1
    A(1,1) = 0.5;
    A(1,2) = 0.5;
    b(1) = 0;
    % BC2
    A(N+2,N+1) = 0.5;
    A(N+2,N+2) = 0.5;
    b(N+2) = CE;
    C = A \ b; % invert matrix
end

function [x, C] = BVP_AD_central_sparse(L, U, Gamma, CE, N)
    dx = L/N;
    x = (-dx/2:dx:L+dx/2)';
    % create diagonal vectors
    r = Gamma / dx^2;
    v1 = (r + 0.5*U/dx) * ones(N+2, 1);
    v2 = -2 * r * ones(N+2, 1);
    v3 = (r -0.5*U/dx) * ones(N+2, 1);
    d = [-1,  0,  1];
    B = [v1, v2, v3];
    % create sparse matrix and enforce BCs
    A = spdiags(B, d, N+2, N+2); 
    A(1,1) = 0.5;     A(1,2) = 0.5;      % BC 1
    A(N+2,N+1) = 0.5; A(N+2,N+2) = 0.5;  % BC 2
    b = zeros(N+2, 1);    
    b(N+2) = CE;
    C = A \ b; % invert matrix
end

function [x, C] = BVP_AD_central_sparse2(L, U, Gamma, CE, N)
    dx = L/N;
    x = (-dx/2:dx:L+dx/2)';
    A = zeros(N+2);
    b = zeros(N+2, 1);
    % interior points
    r1 = 0.5*U/dx;
    r2 = Gamma / dx^2;
    for i = 2:N+1
        A(i, i-1) = r2 + r1;
        A(i, i) = -2 * r2;
        A(i, i+1) = r2 - r1;
        b(i) = 0;
    end
    % BC1
    A(1,1) = 0.5;
    A(1,2) = 0.5;
    b(1) = 0;
    % BC2
    A(N+2,N+1) = 0.5;
    A(N+2,N+2) = 0.5;
    b(N+2) = CE;
    A = sparse(A);
    C = A \ b; % invert matrix
end