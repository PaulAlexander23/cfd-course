import numpy as np
import matplotlib.pyplot as plt


def BVP_AD_central(L, U, Gamma, CE, N):
    dx = L/N
    x = np.linspace(-dx/2, L + dx/2, N + 2)
    A = np.zeros((N + 2, N + 2))
    b = np.zeros(N + 2)
    # Interior points
    r1 = 0.5*U/dx
    r2 = Gamma / dx**2
    for i in range(1, N+1):
        A[i, i-1] = r2 + r1
        A[i, i] = -2 * r2
        A[i, i+1] = r2 - r1
        b[i] = 0
    # BC1
    A[0, 0] = 0.5
    A[0, 1] = 0.5
    b[0] = 0
    # BC2
    A[N+1, N] = 0.5
    A[N+1, N+1] = 0.5
    b[N+1] = CE
    C = np.linalg.solve(A, b)
    return x, C


if __name__ == "__main__":
    L = 1.0
    U = 1.0
    Pe = 10.0
    CE = 1.0
    Gamma = 1.0/Pe

    # exact solution
    def Cex(x):
        return (np.exp(Pe * x/L)-1)/(np.exp(Pe)-1)

    xf = np.linspace(0, L, 1000)

    Ns = [4, 8, 16, 32]
    fig, axarr = plt.subplots(2, 2)
    for m in range(len(Ns)):
        x, Cc = BVP_AD_central(L, U, Gamma, CE, Ns[m])
        dx = x[2] - x[1]
        axarr[m / 2, m % 2].plot(x, Cc, '-^b')
        axarr[m / 2, m % 2].plot(xf, Cex(xf), '-k')
        # axarr[m / 2, m % 2].xlabel('x/L')
        # axarr[m / 2, m % 2].ylabel('C/C_E')
        print 'N = %u; Pe_C= %g' % (Ns[m], U*dx/Gamma)
        axarr[m % 2, m / 2].set_title('N = %u; Pe_C= %g' % (Ns[m], U*dx/Gamma))
    plt.show()
