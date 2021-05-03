import numba
import numpy as np


@numba.njit()
def lda_point(d):   #lda_point(d, SOC):
    if d / (4.0 * np.pi) <= 1.0E-21:
        return 0.0, 0.0

    a0 = (4.0 / 9.0 / np.pi) ** (1.0 / 3.0)
    rs = (3.0 / d) ** (1.0 / 3.0)

    vxp = -2.0 / (np.pi * a0 * rs)
    exxp = 3.0 * vxp / 4.0

    '''
    if SOC:
        a0 = 0.014 / rs
        te = np.sqrt(1.0 + a0 ** 2.0)
        be = np.log(a0 + te)
        vxp *= -0.5 + 1.5 * be / (a0 * te)
        exxp *= 1.0 - 1.5 * ((a0 * te - be) / a0 ** 2.0) ** 2.0
    '''

    if rs > 1.0:
        te = 1.0 + (7.0 / 6.0) * 1.0529 * np.sqrt(rs) + (4.0 / 3.0) * 0.3334 * rs
        be = 1.0 + 1.0529 * np.sqrt(rs) + 0.3334 * rs
        ecp = -0.2846 / be
        vcp = -0.2846 * te / be ** 2.0
    else:
        ecp = 2.0 * ((0.0311 + 0.002 * rs) * np.log(rs) - 0.048 - 0.0116 * rs)
        vcp = 2.0 * ((0.0311 + 2.0 / 3.0 * 0.002 * rs) * np.log(rs) - (0.048 + 0.0311 / 3.0) -
                     (2.0 / 3.0 * 0.0116 + 0.002 / 3.0) * rs)

    return exxp + ecp, vxp + vcp


@numba.njit()
def pbe_point(d, dd):
    dtol = 1.0E-15
    thrd = 1.0 / 3.0

    # Check for zero density and gradient then calculate rs and fs appropriately.
    d1 = max(d, dtol)
    rs = (0.75 / (np.pi * d1)) ** thrd
    fk = (3.0 * np.pi * np.pi * d1) ** thrd
    sk = np.sqrt(4.0 * fk / np.pi)

    if d > dtol:
        s = dd / (d * fk * 2.0)
        t = dd / (d * sk * 2.0)
    else:
        s = 0.0
        t = 0.0

    # Calculate PBE exchange and correlation.
    ex, exd, exdd = pbe_exch(d1, s)
    ec, ecd, ecdd = pbe_corr(d1, rs, t)

    exc = ex + ec
    excd = exd + ecd
    excdd = exdd + ecdd / sk

    return exc, excd, excdd


@numba.njit()
def pbe_exch(rho, s):
    ax = -0.738558766382022405884230032680836
    um = 0.2195149727645171
    uk = 0.8040
    ul = um / uk
    third = 1.0 / 3.0
    ttpi23 = 0.3232409193120513

    # e_x[unif]=ax*rho^(1/3) where ax = -0.75*(3/pi)^(1/3)
    exunif = ax * rho ** third

    # PBE enhancement factor.
    s2 = s * s
    p0 = 1.0 + ul * s2
    p1 = 1.0 / p0

    # fxpbe(s)=1+uk-uk/(1+ul*s*s)
    fxpbe = 1.0 + uk - uk * p1

    # Exchange energy e_x[pbe]=e_x[unif]*fxpbe(s)
    ex = exunif * fxpbe * rho

    # Derivatives of exchange energy.
    fs = 2.0 * um * s * p1 * p1
    exd = 4.0 * exunif * (fxpbe - fs * s) / 3.0
    exdd = ax * 0.5 * fs * ttpi23

    return ex, exd, exdd


@numba.njit()
def pbe_corr(rho, rs, t):
    bet = 0.06672455060314922
    gamma = 0.03109069086965489503494086371273
    delt = bet / gamma
    third = 1.0 / 3.0

    # Get LSD correlation.
    rtrs = np.sqrt(rs)
    ec, eurs = pbe_lsdcor(0.0310907, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294, rtrs)

    # Now do PBE contribution.
    pon = -1.0 * ec / gamma
    q1 = np.exp(pon)
    b = delt / (q1 - 1.0)
    t2 = t * t
    bt2 = b * t2
    q4 = 1.0 + bt2
    q5 = 1.0 + bt2 + bt2 * bt2
    h = gamma * np.log(1.0 + delt * q4 * t2 / q5)
    ect = ec + h

    # PBE correlation energy.
    ec = ect * rho

    # Now do its derivatives.
    q6 = 1.0 / ((gamma * q5 + bet * t2 * q4) * q5)
    q7 = gamma * q6
    ha = -1.0 * bet * t2 * t2 * t2 * b * (1.0 + q4) * q7
    aec = b * b * q1 / bet
    ht = 2.0 * bet * t * (1.0 + 2.0 * bt2) * q7
    h1 = (rs * ha * aec * eurs + 7.0 * t * ht * 0.5) * third
    ecn = third * rs * eurs

    # Set 1st and 2nd derivatives.
    ecd = ect - (ecn + h1)
    ecdd = 0.5 * ht

    return ec, ecd, ecdd


@numba.njit()
def pbe_lsdcor(a, a1, b1, b2, b3, b4, rtrs):
    q0 = -2.0 * a * (1.0 + a1 * rtrs * rtrs)
    q1 = 2.0 * a * rtrs * (b1 + rtrs * (b2 + rtrs * (b3 + b4 * rtrs)))
    q3 = 1.0 + 1.0 / q1
    q2 = np.log(q3)
    gg = q0 * q2
    q4 = 1.0 / rtrs
    q1rs = a * (b1 * q4 + 2.0 * b2 + rtrs * (3.0 * b3 + 4.0 * b4 * rtrs))
    ggrs = -2.0 * a * a1 * q2 - q0 * q1rs / (q1 * (1.0 + q1))

    return gg, ggrs


@numba.njit()
def blyp_point(d, dd):

    ex_b88, exd_b88, exdd_b88 = b88_exchange_point(d, dd)

    ec_lyp, ecd_lyp, ecdd_lyp = lyp_correlation_point(d, dd)

    exc = ex_b88 + ec_lyp
    excd = exd_b88 + ecd_lyp
    excdd = exdd_b88 + ecdd_lyp

    return exc, excd, excdd


@numba.njit()
def b88_exchange_point(d, dd):
    dtol = 1.0E-7

    if d > dtol:
        exc, excd, excdd = sb88_exch(d / 2.0, dd / 2.0)
        exc *= 2.0
    else:
        exc = 0.0
        excd = 0.0
        excdd = 0.0

    return exc, excd, excdd


@numba.njit()
def sb88_exch(d, dd):
    b = 0.0042
    c = 0.930525736349100025
    third = 1.0 / 3.0
    f43 = 4.0 / 3.0

    dd_2      = dd * dd
    d_13      = d ** third
    d_43      = d * d_13
    x         = dd / d_43
    x_2       = x * x
    barsinh_x = b * np.log(x + np.sqrt(1.0 + x_2))
    denom     = 1.0 + 6.0 * x * barsinh_x

    # Exchange energy.
    ex = -1.0 * c * d_43 - b * d_43 * x_2 / denom

    d_73       = d * d_43
    denom_2    = denom * denom
    bdarsinh_x = b / np.sqrt(1.0 + x_2)

    # Derivative of energy w.r.t down density.
    exd = -1.0 * c * f43 * d_13 + f43 * b / d_73 * dd_2 / denom - b / d_43 * dd_2 / denom_2 * 8.0 * (dd / d_73 * barsinh_x + x_2 / d * bdarsinh_x)

    # Derivative of energy w.r.t |grad(down)|
    exdd = -2.0 * b * x / denom + b * dd_2 / d_43 / denom_2 * 6.0 / d_43 * (barsinh_x + x * bdarsinh_x)

    return ex, exd, exdd


@numba.njit()
def lyp_correlation_point(d, dd):
    dtol = 1.0E-20

    if d > dtol:
        ec, ecd, ecdd = lyp_correlation_point_dl(d, dd ** 2.0)
        ex = 0.0
        exd = 0.0
        exdd = 0.0
        exc = ex + ec
        excd = exd + ecd
        excdd = exdd + ecdd * dd / 2.0
    else:
        # Total density is tiny.
        exc = 0.0
        excd = 0.0
        excdd = 0.0

    return exc, excd, excdd


@numba.njit()
def lyp_correlation_point_dl(rhoa1, sigmaaa1):
    tol = 1.0E-20

    rho = max(0.0, rhoa1)

    if rho > tol:
        sigma = max(0.0, sigmaaa1)
        t2 = rho ** (1.0 / 3.0)
        t3 = 1.0 / t2
        t5 = 1.0 + 0.3490 * t3
        t6 = 1.0 / t5
        t9 = 0.2533 * t3
        t10 = np.exp(-1.0 * t9)
        t11 = t10 * t6
        t12 = rho ** 2.0
        t14 = t2 ** 2.0
        t16 = 1.0 / t14 / t12 / rho
        t20 = t3 * t6
        t22 = 0.2611111111111111E1 - 0.9850555555555556E-1 * t3 - 0.13572222222222220 * t20
        t30 = t9 + 0.3490 * t20 - 11.0
        t33 = 0.1148493600075277e2 * t14 * t12 + t22 * sigma +\
              - 0.50 * (0.25E1 - 0.1407222222222222E-1 * t3 - 0.1938888888888889E-1 * t20) * \
              sigma -0.2777777777777778E-1 * t30 * sigma
        t38 = 0.250 * t12 * t33 - 0.45833333333333330 * t12 * sigma
        zk = -0.4918E-1 * t6 * rho - 0.649176E-2 * t11 * t16 * t38
        t45 = t14 * rho
        t49 = t30 / rho * sigma
        t54 = rho * sigma
        t60 = t5 ** 2.0
        t61 = 1.0 / t60
        t64 = t12 ** 2.0
        t66 = 1.0 / t64 / rho
        t81 = 1.0 / t2 / rho
        t83 = t81 * t6
        t85 = 1.0 / t45
        t86 = t85 * t61

        vrhoa = -0.4918E-1 * t6 - 0.649176E-2 * t11 * t16 *\
        (0.50 * rho * t33 + 0.250 * t12 * (0.3062649600200738E2 * t45 - 0.2777777777777778E-1 * t49) -0.250 * t54) +\
        -0.5721273333333333E-2 * t61 * t3 - 0.548120936E-3 * t66 * t10 * t6 * t38 - 0.75520808E-3 * t10 * t61 * t66 * t38 +\
        +0.2380312E-1 * t11 / t14 / t64 * t38 - 0.649176E-2 * t11 * t16 * \
        (0.250 * t12 * ((0.3283518518518519E-1 * t81 + 0.4524074074074074E-1 * t83 - 0.1578901851851852E-1 * t86) * sigma
             -0.50 * (0.4690740740740741E-2 * t81 + 0.6462962962962963E-2 * t83
             -0.2255574074074074E-2 * t86) * sigma - 0.2777777777777778E-1 * (
             -0.8443333333333333E-1 * t81 - 0.11633333333333330 * t83
             +0.4060033333333333E-1 * t86) * sigma + 0.2777777777777778E-1 * t49)
             -0.66666666666666670 * t54)

        vsigmaaa = 0.7213066666666667E-3 * t11 * t85 - 0.2596704E-1 * t11 * t16 * (0.250 * t12 * t22 - 0.66666666666666670 * t12)
    else:
        zk = 0.0
        vrhoa = 0.0
        vsigmaaa = 0.0

    return zk, vrhoa, vsigmaaa