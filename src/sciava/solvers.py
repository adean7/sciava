from numba import njit


def sod_solve(zz, yy, nmin, nmax, h=None):

    # Decide whether integrating from:
    #    left to right ---> isgn = + 1
    # or right to left ---> isgn = - 1
    isgn = (nmax - nmin) // abs(nmax - nmin)

    # Run some test just to be conservative.
    assert isgn == +1 or isgn == -1, 'sod_solve: nmin, nmax out of range'
    if isgn == +1:
        assert nmin > 2, 'sod_solve: nmin, nmax out of range'
        assert nmax <= len(zz), 'sod_solve: nmin, nmax out of range'

    elif isgn == -1:
        assert nmin < (len(zz) - 1), 'sod_solve: nmin, nmax out of range'
        assert nmax >= 1, 'sod_solve: nmin, nmax out of range'

    # Initialise current second derivative of yy (d2y), ycur and yold.
    if h is not None:
        ycur, yold, d2y = sod_solve_der_h(zz, yy, nmin, isgn, h)
    else:
        ycur, yold, d2y = sod_solve_der(zz, yy, nmin, isgn)

    # Begin the integration loop
    if h is not None:
        sod_solve_int_h(zz, yy, nmin, nmax, isgn, h, ycur, yold, d2y)
    else:
        sod_solve_int(zz, yy, nmin, nmax, isgn, ycur, yold, d2y)

@njit()
def sod_solve_der_h(zz, yy, nmin, isgn, h):
    i = nmin - isgn - 1
    d2y = zz[i] * yy[i] + h[i]
    ycur = yy[i] - (zz[i] * yy[i] + h[i]) / 12.0
    i -= isgn
    yold = yy[i] - (zz[i] * yy[i] + h[i]) / 12.0
    return ycur, yold, d2y

@njit()
def sod_solve_der(zz, yy, nmin, isgn):
    i = nmin - isgn - 1
    d2y = zz[i] * yy[i]
    ycur = (1.0 - zz[i] / 12.0) * yy[i]
    i -= isgn
    yold = (1.0 - zz[i] / 12.0) * yy[i]
    return ycur, yold, d2y

@njit()
def sod_solve_int_h(zz, yy, nmin, nmax, isgn, h, ycur, yold, d2y):
    for i in range(nmin - isgn, nmax + isgn, isgn):
        ynew = 2.0 * ycur - yold + d2y
        yy[i] = (ynew + h[i] / 12.0) / (1.0 - zz[i] / 12.0)
        d2y = zz[i] * yy[i] + h[i]
        yold = ycur
        ycur = ynew

@njit()
def sod_solve_int(zz, yy, nmin, nmax, isgn, ycur, yold, d2y):
    for i in range(nmin - isgn, nmax + isgn, isgn):
        ynew = 2.0 * ycur - yold + d2y
        yy[i] = ynew / (1.0 - zz[i] / 12.0)
        d2y = zz[i] * yy[i]
        yold = ycur
        ycur = ynew