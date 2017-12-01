import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    Gamma = 0.005
    L = 1.0
    N = 20
    dx = L/N
    CE = 1.0

    x = np.linspace(-dx/2, L + dx/2, N+2)
    A = np.zeros((N + 2, N + 2))
    b = np.zeros((N + 2, 1))

    r = Gamma/dx**2

    for i in range(1, N+1):
        A[i, i-1] = r
        A[i, i] = - 2 * r
        A[i, i+1] = r
        b[i] = 0

    # BC1
    A[0, 0] = 0.5
    A[0, 1] = 0.5
    b[0] = 0

    # BC2
    A[N + 1, N] = 0.5
    A[N + 1, N + 1] = 0.5
    b[N + 1] = CE

    C = np.linalg.solve(A, b)

    xf = np.linspace(0, L, 1000)
    Cex = xf/L * CE

    plt.plot(x, C, 'ob')
    plt.plot(xf, Cex, '-k')

    plt.xlabel('x/L')
    plt.ylabel('C/CE')
    plt.legend(['numerical', 'exact'])

    plt.show()
