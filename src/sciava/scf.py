import numpy as np

import sciava.algor as algor
import sciava.basis as basis
import sciava.data as data
import sciava.kernels as kernels
import sciava.solvers as solvers



def atomisticCycleSchrod(ab, aeat):
    energy_tol = 1.0E-10
    mixing = 0.4

    # Construct the initial Hamiltonian.
    aeat.init_H(ab)

    # Loop over SCF until converged.
    nit = 0  # Number of iterations.
    max_it = 200  # Maximum number of iterations.
    diff = 0.0
    etoto = 0.0  # Used to detect convergence
    converged = False  # Check if calculation is converged.
    flip = False
    flipo = False

    while not converged and nit < max_it:
        nit += 1

        # Solve the Schrodinger equation.
        aeat.schroedinger_solve(ab)

        # Calculate the charge density, screening, total energy and new potential.
        aeat.set_H(ab, mixing)

        # Determine convergence.
        flip = diff * (aeat.etot - etoto) < 0.0

        diff = aeat.etot - etoto

        if flip and flipo:
            mixing *= 0.9
            flip = False

        converged = abs(diff) < energy_tol
        etoto = aeat.etot
        flipo = flip

    # Calculate the potential energy for this SCF charge density.
    aeat.set_H(ab, mixing=1.0)

    # Output a summary of the calculation.
    print('/===============================================================->')
    print('| Atomic calculation performed for {:>3}: {:<80}'.format(aeat.element, aeat.config))
    print('|---------------------------------------------------------------->')

    # TODO: add in more units
    energy_label = 'eV'
    etot_eV = 27.211386024367247 * aeat.etot  # Convert from Hartrees to eV.
    diff_eV = 27.211386024367247 * diff  # Convert from Hartrees to eV.

    if nit < max_it:
        print('| Converged in {:>3} iterations to an atomic energy of {:>9.3f} {:>3}'.format(nit,
                                                                                             etot_eV,
                                                                                             energy_label))
    else:
        print('| Exceeded max iterations: converged to {:>9.3f} {:>3}'.format(diff_eV, energy_label))
        print('| The atomic energy is {:>9.3f} {:>3}'.format(etot_eV, energy_label))

    print('\===============================================================->')












class AllElectronAtom:
    def __init__(self, currentParams, element: str, Z: int, ab: basis.AllElectronBasis):
        # The elemental symbol.
        self.element = element

        # Allocate the charge.
        self.rhorr = np.zeros(ab.npts)

        # Allocate KED.
        self.kedrr = np.zeros(ab.npts)

        # Allocate the local nuclear potential.
        self.V_nuc = np.zeros(ab.npts)

        # Allocate the local effective potential.
        self.V_eff = np.zeros(ab.npts)
        self.V_tau = np.zeros(ab.npts)

        # Allocate the local valence potential.
        self.V_val = np.zeros(ab.npts)

        # Allocate aeat%norbs.
        self.norbs = np.zeros(ab.lmax + 1, dtype=int)

        # Set the atomic number of this atom.
        self.Z = Z

        # Just so we always have access to it.
        self.cfg_occ = [float(num) for num in data.gs_occ.get(Z, None)]

        # Number of electrons.
        self.nelec = sum(self.cfg_occ)

        for i in range(len(data.gs_occ.get(0))):
            l = data.gs_occ.get(0)[i] % 10
            if self.cfg_occ[i] > 0.0:
                self.norbs[l] = max(self.norbs[l], (data.gs_occ.get(0)[i]-l)/10 - l)
        self.norbmax = np.max(self.norbs)

        # Allocate the real space pseudo orbitals.
        self.orbr = np.zeros((ab.npts, self.norbmax, ab.lmax + 1))

        # Allocate the eigenvalues/energy levels of the orbitals.
        self.eval = np.zeros((self.norbmax, ab.lmax + 1))

        # Allocate the principal quantum number of the orbitals.
        self.orb_n = np.zeros((self.norbmax, ab.lmax + 1), dtype=int)

        # Allocate the occupations of the orbitals.
        self.occ = np.zeros((self.norbmax, ab.lmax + 1), dtype=int)

        # Total energy.
        self.etot = 0.0

        # Finite nucleus?
        self.fin_nuc = True if currentParams.atomicSolver in ['SH', 'KH'] else False

        # Configuration of electrons. 1s2, 2s1, etc.
        self.config = ' '

        # Get the DFT theory.
        self.theory = currentParams.theory

    def init_H(self, ab: basis.AllElectronBasis):
        # Set orbital occupancies.
        self.set_occ(ab)

        # Establish an analytic starting potential.
        for n in range(ab.npts):
            arg = ab.r[n] / ab.r[29]
            if arg < 20.0:
                tt = np.exp(-arg)
            else:
                tt = 0.0

            # Use a finite nucleus?
            if self.fin_nuc:
                self.V_nuc[n] = -self.Z * (1.0 - tt)
            else:
                self.V_nuc[n] = -self.Z # Use a point nucleus.

        drng = 0.675
        zeff = min(float(self.Z), self.nelec) - 0.5
        yy   = zeff ** 0.4
        for n in range(ab.npts):
            tt = min(100.0, ab.r[n]/drng)
            tt = zeff * (1.0 - 1.0 / (drng * yy * (np.exp(tt) - 1.0) + 1.0))
            self.V_eff[n] = tt + self.V_nuc[n]

    def set_occ(self, ab: basis.AllElectronBasis):
        if self.Z == 0:
            # self.occ already set to ' '
            return

        # Analyse the configuration.
        self.norbs = np.zeros(ab.lmax + 1, dtype=int)
        self.config = data.gs_config.get(self.Z)
        no2cf = np.zeros((10, ab.lmax + 1), dtype=int)
        for i in range(len(self.cfg_occ)):
            if self.cfg_occ[i] > 0.0:
                l = data.gs_occ.get(0)[i] % 10
                self.norbs[l] += 1
                no2cf[self.norbs[l]-1,l] = i+1

        # Set the occupations.
        num = [0, 0, 0, 0] # Number of s, p, d, f angular momentum channels.
        for i in range(len(self.cfg_occ)):
            if self.cfg_occ[i] > 0:
                l = data.gs_occ.get(0)[i] % 10
                num[l] += 1
                assert num[l] <= self.norbmax, 'num(l) > aeat.norbmx'
                self.occ[num[l]-1,l] = self.cfg_occ[i]

        # Hydrogenic(ish) energy levels as a first guess (a fudged screening).
        nmax = 0
        for l in range(ab.lmax+1):
            for nob in range(self.norbs[l]):
                n = data.gs_occ.get(0)[no2cf[nob,l]-1] // 10
                if n > nmax:
                    nmax = n

        for l in range(ab.lmax+1):
            for nob in range(self.norbs[l]):
                n = data.gs_occ.get(0)[no2cf[nob,l]-1] // 10
                self.orb_n[nob,l] = n
                self.eval[nob,l]  = -(1.0 + self.Z - 0.9 * self.Z / float(nmax) * float(n)) ** 2.0 / float(n) ** 2.0 / 2.0

    def schroedinger_solve(self, ab: basis.AllElectronBasis):
        #      integer       :: i,l,no,niter
        #      integer       :: nodes,nodes_target
        #      integer       :: nctp,ninf
        #      integer       :: max_solve
        #      real(kind=dp) :: solver_thresh
        #      real(kind=dp) :: ecur,emin,emax
        #      real(kind=dp) :: vzero,d2p,disc
        #      real(kind=dp) :: factor,norm,decur,yyctp
        #      real(kind=dp) :: zz(ab%npts),yy(ab%npts),vme(ab%npts)
        #      real(kind=dp) :: yval(-1:+1)
        #      real(kind=dp) :: aa,bb,Z

        assert self.fin_nuc, 'schroedinger_solve: Must be a finite nucleus'

        solver_thresh = 1.0E-10
        max_solve = 1000
        yval = np.zeros(3)

        for l in range(ab.lmax+1):

            for no in range(self.norbs[l]):

                ecur = self.eval[no,l] * 2.0              # Set the current eigenvalue (convert to Ry)
                nodes_target = self.orb_n[no,l] - l - 1   # Set the expected number of nodes

                emin = -1.0E10  # Set initial upper and
                emax =  1.0     #   lower bounds for the eigenvalue

                niter = 1                  # Reset the iteration counter
                solver_converged = False

                yy  = np.zeros(ab.npts)      # Zero the wavefunction
                zz  = np.zeros(ab.npts)
                vme = np.zeros(ab.npts)

                # ------------------------------------------------------
                # I t e r a t i v e  s o l u t i o n  f o r  e n e r g y
                # ------------------------------------------------------

                while not solver_converged and niter < max_solve:
                    niter += 1

                    # Set ( v(r) - e ) Note: convert V_eff to Ry
                    vme[1:] = 2.0 * self.V_eff[1:] / ab.r[1:] - ecur

                    # Compute the potential for this angular momentum and energy.
                    zz[0] = 0.0
                    zz[1:] = ab.b ** 2.0 / 4.0 + ab.rab[1:] ** 2.0 * ((2.0 * self.V_eff[1:] / ab.r[1:]) +
                              (float(l*(l+1)) / ab.r[1:] ** 2.0) - ecur)

                    # Find the classical turning point...
                    for nctp in range(ab.npts-1, 8, -1):
                        if vme[nctp] < 0.0:
                            break

                    assert nctp != 8, 'schroedinger_solve: no classical turning point found'
                    assert nctp <= ab.npts-11, 'schroedinger_solve: classical turning point too close to ab.npts'

                    # ... and practical infinity
                    for ninf in range(nctp + 10, ab.npts):
                        if (vme[ninf] * (ab.r[ninf] - ab.r[nctp]) ** 2.0) > (np.log(solver_thresh) ** 2.0):
                            break

                    assert ninf <= ab.npts, 'schroedinger_solve: ninf too big for this mesh'

                    # ============================================
                    # =             O U T W A R D S              =
                    # ============================================

                    # S o l v e  t h e  h o m o g e n e o u s  e q u a t i o n

                    # Compute the starting value of yy according to
                    # yy(r) = r ** (l+1) * ( 1 + aa * r + bb * r**2 + ... )
                    #   where:
                    # aa = - z / (l+1)
                    # bb = ( 2 * z * z / (l+1) + v(0) - e ) / (4*l + 6)

                    '''
                    Z = 0.0

                    vzero = 2.0 * self.V_eff[1] / ab.r[1]
                    aa = - Z / float(l + 1)
                    bb = ((2.0 * Z ** 2.0 / float(l+1)) + vzero - ecur) / float(4 * l + 6)

                    yy[0:3] = ab.r[0:3] ** float(l + 1) * (1.0 + aa * ab.r[0:3] + bb * ab.r[0:3] ** 2.0)
                    '''

                    vzero = 2.0 * self.V_eff[1] / ab.r[1]
                    bb = (vzero - ecur) / float(4 * l + 6)

                    yy[0:3] = ab.r[0:3] ** float(l + 1) * (1.0 + bb * ab.r[0:3] ** 2.0)

                    # Convert to the phi function.
                    yy[0] = yy[0] / ab.sqr[0]
                    yy[1] = yy[1] / ab.sqr[1]
                    yy[2] = yy[2] / ab.sqr[2]

                    # Integrate the schroedinger equation outwards.
                    solvers.sod_solve(zz, yy, 4, nctp)

                    # Save the value at nctp.
                    yyctp = yy[nctp]

                    # ============================================
                    # =              I N W A R D S               =
                    # ============================================

                    # Start up of wavefunction at practical infinity.
                    yy[ninf]   = np.exp(-np.sqrt(zz[ninf]  ) * float(ninf-nctp  ))
                    yy[ninf-1] = np.exp(-np.sqrt(zz[ninf-1]) * float(ninf-nctp-1))

                    # I n w a r d   i n t e g r a t i o n   t o   n c t p
                    solvers.sod_solve(zz, yy, ninf-2, nctp)

                    # ============================================
                    # =     M A T C H  A N D  O P T I M I S E    =
                    # ============================================

                    # Rescale tail to make continuous at nctp.
                    factor = yyctp / yy[nctp]
                    yy[nctp:ninf] = yy[nctp:ninf] * factor

                    # Check that the number of nodes is OK
                    nodes = algor.number_of_nodes(yy, 1, ninf)

                    if nodes < nodes_target: # Energy is too low
                        emin = ecur
                        if ecur * 0.97 > emax:
                            ecur  = 0.5 * (ecur + emax)
                        else:
                            ecur *= 0.97
                        continue

                    if nodes > nodes_target: # Energy is too high
                        emax = ecur
                        if ecur * 1.03 < emin:
                            ecur  = 0.5 * (ecur + emin)
                        else:
                            ecur *= 1.03
                        continue

                    # Find normalisation of wavefunction using Simpson's rule.
                    norm = basis.radin(ab.npts, ab.rab, (yy * ab.sqr) ** 2.0)

                    # Compute the Noumerov discontinuity in the yy function.
                    for i in range(len(yval)):
                        yval[i] = yy[nctp+i-1] * (1.0 - zz[nctp+i-1] / 12.0)
                    d2p  = zz[nctp] * yyctp
                    disc = -yval[0] + 2.0 * yval[1] - yval[2] + d2p

                    # Variational estimate for the change in ecur.
                    decur = yyctp * disc * ab.sqr[nctp] ** 2.0 / (norm * ab.rab[nctp])

                    # To prevent convergence problems:
                    # do not allow decur to exceed 20% of |ecur|
                    # do not allow decur to exceed 70% of distance to emin or emax
                    if decur > 0.0:
                        emin = ecur
                        decur = min(decur, -0.2 * ecur, 0.7 * (emax - ecur))
                    else:
                        emax = ecur
                        decur = -min(-decur, -0.2 * ecur, 0.7 * (ecur - emin))

                    # Test to see whether eigenvalue converged.
                    solver_converged = abs(decur) < solver_thresh

                    ecur += decur

                    # Check that the iterative loop is not about to terminate.
                    # Eigenfunction has not converged in allowed number of iterations.
                    assert niter < max_solve, 'schrodinger_solve: could not converge'

                # Update the energy.
                self.eval[no,l] = ecur / 2.0    # Note: convert to Hartree

                # Update the wavefunction array.
                self.orbr[:,no,l] = yy * ab.sqr / np.sqrt(norm)

    def set_H(self, ab: basis.AllElectronBasis, mixing: float):
        # Calculate the charge density.

        self.V_val = np.zeros(ab.npts)
        self.rhorr = np.zeros(ab.npts)

        for l in range(ab.lmax+1):
            for no in range(self.norbs[l]):
                self.rhorr += self.occ[no,l] * self.orbr[:,no,l] ** 2.0

        work  = np.zeros(ab.npts)
        #work2 = np.zeros(ab.npts)

        '''
        work[1:] = self.rhorr[1:] / ab.r[1:]

        basis.r_second_derivative(ab.npts, ab.rab, ab.r, f=work, df_dr=work2, d2f_dr2=self.kedrr)
        self.kedrr *= ab.r

        for l in range(ab.lmax+1):
            for no in range(self.norbs[l]):
                basis.r_second_derivative(ab.npts, ab.rab, ab.r, f=self.orbr[:,no,l], df_dr=work, d2f_dr2=work2)

                work[1:] = work2[1:] - l * (l+1) * self.orbr[1:,no,l] / ab.r[1:] ** 2.0

                self.kedrr -= self.occ[no,l] * 2.0 * work * self.orbr[:,no,l]

        self.kedrr *= 0.25
        '''

        # No boundary conditions in standalone so no need to put rho_sps onto rhorr_val.

        # Calculate the corresponding effective local potential.
        V_eff0 = self.V_eff      # Save the old potential for mixing.
        V_tau0 = self.V_tau      # Save the old potential for mixing.

        self.V_eff = self.V_nuc.copy()  # The external local potential xr

        if any(self.rhorr > 0.0):
            Ehar = self.hartree(ab, work)

            self.V_eff += work * ab.r # The Hartree potential.

        # No boundary conditions in standalone so no need to deal with self.rhorr_val

        de_dk = np.zeros(ab.npts)
        Exc = self.xc(ab, work, de_dk=de_dk)

        self.V_eff += work * ab.r   # The XC potential.
        self.V_tau = de_dk

        # Mix the old and new potentials.
        self.V_eff = mixing * self.V_eff + (1.0 - mixing) * V_eff0
        self.V_tau = mixing * self.V_tau + (1.0 - mixing) * V_tau0

        # Calculate the total energy.
        self.etot = 0.0
        for l in range(ab.lmax+1):
            for no in range(self.norbs[l]):
                self.etot += self.occ[no,l] * self.eval[no,l]

        self.etot += Ehar + Exc

    def hartree(self, ab, vhar, nrc=None, Vrc=None):
        work1 = np.zeros(ab.npts)
        work2 = np.zeros(ab.npts)
        func = np.zeros(ab.npts)

        v1 = basis.radin(ab.npts, ab.rab, 2.0 * self.rhorr, inter=work1)

        func[0] = 0.0
        func[1:] = 2.0 * self.rhorr[1:] / ab.r[1:]

        v2 = basis.radin(ab.npts, ab.rab, func, inter=work2)

        vhar[0] = 0.0
        vhar[1:] = work1[1:] / ab.r[1:] - work2[1:] + v2

        ehar = -1.0 * basis.radin(ab.npts, ab.rab, vhar * self.rhorr) / 2.0

        vhar /= 2.0  # Convert to Hartrees.
        ehar /= 2.0  # Convert to Hartrees.

        if nrc is not None:
            assert Vrc is not None, 'aeat_hartree: Vrc must be supplied if nrc is not None'

            work1[0] = 0.0
            work1[1:] = -1.0 * self.rhorr[1:] / ab.r[1:]

            v1 = basis.radin(ab.npts, ab.rab, work1, inter=func)
            v1 = basis.radin(ab.npts, ab.rab, func, inter=work1)

            vhar[0] = 0.0
            vhar[1:] = work1[1:] / ab.r[1:] + Vrc / ab.r[nrc] - work1[nrc] / ab.r[nrc]

            ehar = -1.0 * basis.radin(nrc, ab.rab, vhar * self.rhorr) / 2.0

        return ehar

    def xc(self, ab, vxc, kedrr=None, de_dk=None):

        # Check sizes are sufficient.
        assert len(self.rhorr) >= ab.npts, 'Error, rhorr must have at least ab.npts in ae_xc'
        assert len(vxc) >= ab.npts, 'Error, vxc must have at least ab.npts in ae_xc'

        # Set rho from rhorr.
        rho = np.zeros(ab.npts)
        rho[1:] = self.rhorr[1:] / ab.r[1:] ** 2.0
        rho[0] = rho[2] - (rho[3] - rho[1]) / 2.0 / ab.rab[2] * (ab.r[2] - ab.r[0])         # Extrapolate back.

        if kedrr is not None:
            assert len(kedrr) >= ab.npts, 'Error, kedrr must have at least ab.npts in ae_xc'
            ked = np.zeros(ab.npts)
            ked[1:] = kedrr[1:] / ab.r[1:] ** 2.0
            ked[0] = ked[2] - (ked[3] - ked[2]) / 2.0 / ab.rab[2] * (ab.r[2] - ab.r[0])     # Extrapolate back.

        if de_dk is not None:
            de_dk = np.zeros(ab.npts)

        # Fit an exponential tail.
        rho_thresh = 1.0E-10

        iexp = ab.npts
        for n in range(ab.npts):
            if rho[n] > rho_thresh:
                iexp = n
        iexp = min(iexp, ab.npts - 1)

        grad = np.zeros(ab.npts)
        basis.real_derivative(ab.npts, ab.rab, rho, grad)

        bexp = (grad[iexp] / rho[iexp] + grad[iexp-1] / rho[iexp-1] + grad[iexp+1] / rho[iexp+1]) / 3.0
        aexp = rho[iexp] / np.exp(bexp * ab.r[iexp])

        rho[iexp:] = aexp * np.exp(bexp * ab.r[iexp:])

        work = np.zeros(ab.npts)
        dv_dn = np.zeros(ab.npts)

        if self.theory == 'LDA':
            # The LDA of Ceperley-Alder

            for n in range(ab.npts):
                work[n], vxc[n] = kernels.lda_point(rho[n])

            # Evaluate the XC energy.
            exc = -1.0 * basis.radin(ab.npts, ab.rab, (vxc - work) * self.rhorr)

            vxc /= 2.0  # Convert to Hartrees.
            exc /= 2.0  # Convert to Hartrees.

        elif self.theory == 'PBE':
            rho /= 4.0 * np.pi      # The XC routines derive from a non-radial form

            basis.real_derivative(ab.npts, ab.rab, rho, grad)

            mod_grad = np.abs(grad)

            for n in range(ab.npts):
                work[n], vxc[n], dv_dn[n] = kernels.pbe_point(rho[n], mod_grad[n])

                if mod_grad[n] > 0.0:
                    dv_dn[n] *= -1.0 * grad[n] / mod_grad[n]
                else:
                    dv_dn[n] = 0.0

            # Add the gradient corrections to the potential.
            basis.real_derivative(ab.npts, ab.rab, dv_dn, grad)

            grad[1:] += 2.0 * dv_dn[1:] / ab.r[1:]

            vxc += grad

            # Calculate the eigenvalue correction energy.
            exc = -1.0 * basis.radin(ab.npts, ab.rab, vxc * self.rhorr - work * ab.r ** 2.0 * 4.0 * np.pi)

        elif self.theory == 'BLYP':
            rho /= 4.0 * np.pi  # The XC routines derive from a non-radial form

            basis.real_derivative(ab.npts, ab.rab, rho, grad)

            mod_grad = np.abs(grad)

            for n in range(ab.npts):
                work[n], vxc[n], dv_dn[n] = kernels.blyp_point(rho[n], mod_grad[n])

                if mod_grad[n] > 0.0:
                    dv_dn[n] *= -1.0 * grad[n] / mod_grad[n]
                else:
                    dv_dn[n] = 0.0

            # Add the gradient corrections to the potential.
            basis.real_derivative(ab.npts, ab.rab, dv_dn, grad)

            grad[1:] += 2.0 * dv_dn[1:] / ab.r[1:]

            vxc += grad

            # Calculate the eigenvalue correction energy.
            exc = -1.0 * basis.radin(ab.npts, ab.rab, vxc * self.rhorr - work * ab.r ** 2.0 * 4.0 * np.pi)

        else:
            assert False, 'XC functional theory {} not known/implemented.'.format(self.theory)

        return exc
