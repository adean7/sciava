import numba
import numpy as np


@numba.njit()
def real_derivative(npts, rab, f, df_dr):
    twothird = 2.0 / 3.0

    for i in range(2, npts-2):
        f1 = (f[i+1] - f[i-1]) * twothird
        f2 = (f[i+2] - f[i-2]) / 12.0
        df_dr[i] = (f1 - f2)

    df_dr[1] = 2.0 * df_dr[2] - df_dr[3]
    df_dr[0] = 0.0 # By symmetry.
    df_dr[-2] = 2.0 * df_dr[-3] - df_dr[-4]
    df_dr[-1] = 2.0 * df_dr[-2] - df_dr[-3]

    # Now multiply di/dr.
    for i in range(npts):
        df_dr[i] /= rab[i]


@numba.njit()
def r_second_derivative(npts, rab, r, f, df_dr, d2f_dr2):
    a = [-12.0, 108.0, -540.0,    0.0, 540.0, -108.0, 12.0]
    b = [  4.0, -54.0,  540.0, -980.0, 540.0,  -54.0,  4.0]

    aa = r[npts-1] * rab[0] / (rab[npts-1] - rab[0])

    for i in range(3, npts-3):
        # The first derivative on the log grid
        sum1 = 0.0
        for j in range(7):
            sum1 += a[j]*f[j-3+i]
        sum1 /= 720.0 * rab[i]

        df_dr[i] = sum1

        # The second derivative on the reg grid.
        sum2 = 0.0
        for j in range(7):
            sum2 += b[j] * f[j-3+i]
        sum2 /= 360.0

        d2f_dr2[i] = sum2 / rab[i] ** 2.0 - sum1 / (r[i] + aa)

    df_dr[-3:] = 0.0
    df_dr[:3] = r[:3] * df_dr[3] / r[3]

    d2f_dr2[-3:] = 0.0
    d2f_dr2[:3] = r[:3] * d2f_dr2[3] / r[3]


@numba.njit()
def radin(npts, rab, func, inter=None):
    # If an even number of logarithmic grid points are supplied, use all but the last.
    npts = npts if npts % 2 == 1 else npts-1

    integral = 0.0
    f3 = func[0] * rab[0]

    if inter is None:
        for i in range(1, npts, 2):
            f1 = f3
            f2 = func[i] * rab[i]
            f3 = func[i+1] * rab[i+1]
            integral += (f1 + 4.0 * f2 + f3) / 3.0

    else:
        inter[0] = 0.0
        r12 = 1.0 / 12.0
        for i in range(1, npts, 2):
            f1 = f3
            f2 = func[i] * rab[i] * r12
            f3 = func[i+1] * rab[i+1] * r12
            integral +=  5.0 * f1 + 8.0 * f2 - f3
            inter[i] = integral
            integral += -1.0 * f1 + 8.0 * f2 + 5.0 * f3
            inter[i+1] = integral

    return integral



class AllElectronBasis:
    npts : int      = None # Number of points in the logarithmic grid
    rmax : float    = None # The maximum extent of the logarithmic grid
    a    : float    = None # A parameter defining the logarithmic grid
    b    : float    = None # A parameter defining the logarithmic grid
    r    : np.array = None # The logarithmic grid
    rab  : np.array = None # Logarithmic grid array (for integration)
    sqr  : np.array = None # Logarithmic grid array
    rwgt : np.array = None # Logarithmic grid weight array (for integration)
    lmax : int      = None # The maximum angular momentum channel

    def __init__(self, Z):
        """
        Initialise the basis for the all electron atomic calculations
        :param Z:  The atomic number of the species
        """

        aasf = 6.0
        bbsf = 240.0

        self.a = np.exp(-aasf) / float(Z)
        self.b = 1.0 / bbsf

        self.rmax = 400.0  # Set the extent of the logarithmic grid

        self.npts = 2 + int(np.log(self.rmax / self.a + 1.0) / self.b)  # Set the size of the grid
        self.npts = self.npts + 1 if self.npts % 2 == 0 else self.npts  # Convert to an odd number for Simpson's rule

        # Initialise arrays.
        self.r = np.zeros(self.npts)
        self.rab = np.zeros(self.npts)
        self.sqr = np.zeros(self.npts)
        self.rwgt = np.zeros(self.npts)

        # Set the grid arrays
        for i in range(self.npts):
            self.r[i] = self.a * (np.exp(self.b * float(i)) - 1.0)
            self.rab[i] = self.b * (self.r[i] + self.a)
            self.sqr[i] = np.exp(self.b * float(i + 1) / 2.0)

        # Set the integration weights
        self.rwgt[0] = 1.0 / 3.0  # Initial point.
        for i in range(1, self.npts - 1, 2):  # Even points.
            self.rwgt[i] = self.rab[i] * 4.0 / 3.0

        for i in range(2, self.npts - 2, 2):  # Odd points (excluding the final point).
            self.rwgt[i] = self.rab[i] * 2.0 / 3.0
        self.rwgt[self.npts - 1] = 1.0 / 3.0  # Final point.

        # Set the maximum angular momentum
        self.lmax = 3
