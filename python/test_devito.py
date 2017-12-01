import time
from argparse import ArgumentParser

import numpy as np
import sympy

from devito import Grid, Eq, Operator, TimeFunction, Function
from devito.logger import log


def main():

    Gamma = 0.005
    U = 0.05
    L = 1
    CE = 1

    nx = 100
    x = np.linspace(0, L, nx)

    grid = Grid((nx, ))
    print grid

    C = Function(name='C', grid=grid, space_order=2)
    print C

    eqn = Eq(Gamma * C.dx2 - U * C.dx)
    print eqn

    stencil = sympy.solve(eqn, C)[0]
    print stencil

    bc = Function(name='bc', grid=grid, space_order=2)
    bc.data[0] = 0
    bc.data[nx] = CE
    bc
    print bc

    op = Operator(expressions=[eqn] + bc)
    print op

    # Execute the generated Devito stencil operator
    # tstart = time.time()
    # op.apply(u=u)
    # runtime = time.time() - tstart


if __name__ == "__main__":
    main()
