#
# This class is imported from the abipy package: https://github.com/abinit/abipy
# 
#


# coding: utf-8
"""
Shankland-Koelling-Wood Fourier interpolation scheme.
For the theoretical background see :cite:`Euwema1969,Koelling1986,Pickett1988,Madsen2006`.
"""
import itertools
import numpy as np
import scipy
import time
from collections import deque
from monty.termcolor import cprint
from monty.collections import dict2namedtuple

class SkwInterpolator():
    """
    This object implements the Shankland-Koelling-Wood Fourier interpolation scheme.
    It can be used to interpolate functions in k-space with the periodicity of the
    reciprocal lattice and satisfying F(k) = F(Sk) for each rotation S
    belonging to the point group of the crystal. For readability reason,
    the names of the variables are chosen assuming we are interpolating electronic eigenvalues
    but the same object can be used to interpolate other quantities. Just set the first dimension to 1.
    """

    def __init__(self, lpratio, kpts, eigens, fermie, nelect, cell, symrel, has_timrev,
                 filter_params=None, verbose=1):
        """
        Args:
            lpratio: Ratio between the number of star-functions and the number of ab-initio k-points.
                5-10 should be OK in many systems, larger values may be required for accurate derivatives.
            kpts: numpy array with the [nkpt, 3] ab-initio k-points in reduced coordinates.
            eigens: numpy array with the ab-initio energies. shape [nsppol, nkpt, nband].
            fermie: Fermi energy in eV.
            nelect: Number of electrons in the unit cell
            cell: (lattice, positions, numbers)
                lattice: numpy array with direct lattice vectors along the rows.
                positions: atomic positions in reduced coordinates.
                numbers: Atomic number for each atom in positions.
            symrel: [nsym, 3, 3] numpy array with the (ferromagnetic) symmetry operations of the direct lattice
                in reduced coordinates. anti-ferromagnetic symmetries (if any) should be removed by the caller.
            has_timrev: True is time-reversal can be used.
            filter_params: List with parameters used to filter high-frequency components (Eq 9 of PhysRevB.61.1639)
                First item gives rcut, second item sigma. Ignored if None.
            verbose: Verbosity level.
        """
        self.verbose = verbose
        self.cell = cell
        lattice = self.cell[0]
        self.original_fermie = fermie
        self.interpolated_fermie = self.original_fermie
        self.nelect = nelect
        self.has_timrev = has_timrev

        # iscomplexobj is used to handle lifetimes.
        eigens = np.atleast_3d(eigens)
        self.iscomplexobj = np.iscomplexobj(eigens)
        self.nsppol, self.nkpt, self.nband = eigens.shape

        if len(kpts) != self.nkpt:
            raise ValueError("Second dimension of eigens should be %d but got array of shape: %s" %
                (len(kpts), eigens.shape))
        if self.nkpt == 1:
            raise ValueError("Interpolation algorithm requires nkpt > 1")

        rprimd = np.asarray(lattice).T
        self.rmet = np.matmul(rprimd.T, rprimd)

        # Find point group operations.
        symrel = np.reshape(symrel, (-1, 3, 3))
        self.ptg_symrel, self.ptg_symrec, has_inversion = extract_point_group(symrel, has_timrev)
        self.ptg_nsym = len(self.ptg_symrel)
        if self.verbose:
            print("Found", self.ptg_nsym, "symmetries in point group")

        # Find nrwant star points.
        self.lpratio = lpratio = int(lpratio)
        if lpratio <= 1:
            raise ValueError("lpratio must be > 1 but got %s" % lpratio)

        nrwant = lpratio * self.nkpt
        fact = 1/2 if has_inversion else 1
        rmax = int((1.0 + (lpratio * self.nkpt * self.ptg_nsym * fact) / 2.0) ** (1/3.)) * np.ones(3, dtype=int)
        #rmax = int((1.0 + (lpratio * self.nkpt) / 2.0) ** (1/3.)) * np.ones(3, dtype=int)

        while True:
            self.rpts, r2vals, ok = self._find_rstar_gen(nrwant, rmax)
            self.nr = len(self.rpts)
            if ok:
                break
            else:
                print("rmax: ", rmax," was not large enough to find", nrwant, "R-star points.")
                rmax *= 2
                print("Will try again with enlarged rmax:", rmax)

        print("Using:", self.nr, "star-functions. nstars/nk:", self.nr / self.nkpt)

        # Compute (inverse) roughness function.
        c1, c2 = 0.25, 0.25
        r2min = r2vals[1]
        inv_rhor = 1.0 / ((1.0 - c1 * r2vals / r2min) ** 2 + c2 * (r2vals / r2min) ** 3)

        # Construct star functions for the ab-initio k-points.
        nsppol, nband, nkpt, nr = self.nsppol, self.nband, self.nkpt, self.nr
        self.skr = np.empty((nkpt, nr), dtype=complex)
        for ik, kpt in enumerate(kpts):
            self.skr[ik] = self.get_stark(kpt)

        # Build H(k,k') matrix (Hermitian)
        hmat = np.empty((nkpt-1, nkpt-1), dtype=complex)
        for jk in range(nkpt-1):
            v_jkr = self.skr[jk, 1:] - self.skr[nkpt-1, 1:]
            #for ik in range(jk + 1):
            for ik in range(nkpt-1):
                v_ikr = inv_rhor[1:] * (self.skr[ik, 1:] - self.skr[nkpt-1, 1:])
                hmat[ik, jk] = np.vdot(v_jkr, v_ikr)
                if ik == jk: hmat[ik, jk] = hmat[ik, jk].real

        # Solving system of linear equations to get lambda coeffients (eq. 10 of PRB 38 2721)..."
        de_kbs = np.empty((nkpt-1, nband, nsppol), dtype=complex)
        for spin in range(nsppol):
            for ib in range(nband):
                de_kbs[:, ib, spin] = eigens[spin, 0:nkpt-1, ib] - eigens[spin, nkpt-1, ib]

        # Solve all bands and spins at once
        # FIXME: Portability problem with scipy 0.19 in which linalg.solve wraps the expert drivers
        # http://scipy.github.io/devdocs/release.0.19.0.html#foreign-function-interface-improvements
        if scipy.__version__ == "0.19.0":
            import warnings
            warnings.warn("linalg.solve in scipy 0.19.0 gives weird results. Use at your own risk!!!")

        try:
            lmb_kbs = scipy.linalg.solve(hmat, np.reshape(de_kbs, (-1, nband * nsppol)))
            #lmb_kbs = scipy.linalg.solve(hmat, np.reshape(de_kbs, (-1, nband * nsppol)),
                    #sym_pos=True, lower=False, overwrite_a=True, overwrite_b=True, check_finite=False)
                    #sym_pos=False, lower=False, overwrite_a=False, overwrite_b=False, check_finite=True)

        except scipy.linalg.LinAlgError as exc:
            print("Cannot solve system of linear equations to get lambda coeffients (eq. 10 of PRB 38 2721)")
            print("This usually happens when there are symmetrical k-points passed to the interpolator.")
            raise exc

        lmb_kbs = np.reshape(lmb_kbs, (-1, nband, nsppol))

        # Compute coefficients.
        self.coefs = np.empty((nsppol, nband, nr), dtype=complex)
        for spin in range(nsppol):
            for ib in range(nband):
                for ir in range(1, nr):
                    self.coefs[spin, ib, ir] = inv_rhor[ir] \
                        * np.vdot(self.skr[:nkpt-1, ir] - self.skr[nkpt-1, ir], lmb_kbs[:nkpt-1, ib, spin])

                self.coefs[spin, ib, 0] = eigens[spin, nkpt-1, ib] \
                    - np.dot(self.coefs[spin, ib, 1:nr], self.skr[nkpt-1, 1:nr])

        # Filter high-frequency.
        self.rcut, self.rsigma = None, None
        if filter_params is not None:
            self.rcut = filter_params[0] * np.sqrt(r2vals[-1])
            self.rsigma = rsigma = filter_params[1]
            if self.verbose:
                print("Applying filter (Eq 9 of PhysRevB.61.1639) with rcut:", self.rcut, ", rsigma", self.rsigma)
            from scipy.special import erfc
            for ir in range(1, nr):
                self.coefs[:, :, ir] *= 0.5 * erfc((np.sqrt(r2vals[ir]) - self.rcut) / self.rsigma)

        # Prepare workspace arrays for star functions.
        self.cached_kpt = np.ones(3) * np.inf
        self.cached_kpt_dk1 = np.ones(3) * np.inf
        self.cached_kpt_dk2 = np.ones(3) * np.inf

        # Compare ab-initio data with interpolated results.
        mae = 0.0
        for spin in range(nsppol):
            for ik, kpt in enumerate(kpts):
                skw_eb = self.eval_sk(spin, kpt)
                mae += np.abs(eigens[spin, ik] - skw_eb).sum()
                if self.verbose >= 10:
                    # print interpolated eigenvales
                    for band in range(self.nband):
                        e0 = eigens[spin, ik, band]
                        eskw = skw_eb[band]
                        print("spin", spin, "band", band, "ikpt", ik, "e0", e0, "eskw", eskw, "diff", e0 - eskw)

        mae *= 1e3 / (nsppol * nkpt * nband)
        if np.isnan(mae) or np.isinf(mae) or mae > 1000:
            raise RuntimeError("Interpolation went bananas! mae = %s" % mae)

        warn = mae > 10.0
        cprint("FIT vs input data: Mean Absolute Error= %.3e (meV)" % mae, color="red" if warn else "green")
        if warn:
            # Issue warning if error too large.
            cprint("Large error in SKW interpolation!", "red")
            cprint("MAE:", mae, "[meV]", "red")

        self.mae = mae

    def __str__(self):
        return self.to_string()

    def to_string(self, verbose: int = 0) -> str:
        """String representation."""
        lines = []
        app = lines.append
        app("nsppol: %s, nband: %s" % (self.nsppol, self.nband))
        app("Number of operations in point-group: %s with time-reversal: %s" % (self.ptg_nsym, self.has_timrev))
        app("Number of ab-initio k-points: %s" % self.nkpt)
        app("Number of star functions: %s [nstars/nk = %s]" % (self.nr, (self.nr / self.nkpt)))
        if self.rcut is not None:
            app("Fourier filter (Eq 9 of PhysRevB.61.1639) with rcut: %s, rsigma: %s" % (self.rcut, self.rsigma))
        app("Comparison between ab-initio data and fit gave Mean Absolute Error: %s [meV]" % self.mae)

        return "\n".join(lines)

    def interp_kpts(self, kfrac_coords, dk1=False, dk2=False):
        """
        Interpolate energies on an arbitrary set of k-points. Optionally, compute
        gradients and Hessian matrices.

        Args:
            kfrac_coords: K-points in reduced coordinates.
            dk1 (bool): True if gradient is wanted.
            dk2 (bool): True to compute 2nd order derivatives.

        Return:
            namedtuple with:
            interpolated energies in eigens[nsppol, len(kfrac_coords), nband]
            gradient in dedk[self.nsppol, len(kfrac_coords), self.nband, 3))
            hessian in dedk2[self.nsppol, len(kfrac_coords), self.nband, 3, 3))

            gradient and hessian are set to None if not computed.
        """
        start = time.time()

        kfrac_coords = np.reshape(kfrac_coords, (-1, 3))
        new_nkpt = len(kfrac_coords)
        new_eigens = np.empty((self.nsppol, new_nkpt, self.nband))

        dedk = None if not dk1 else np.empty((self.nsppol, new_nkpt, self.nband, 3))
        dedk2 = None if not dk2 else np.empty((self.nsppol, new_nkpt, self.nband, 3, 3))

        der1, der2 = None, None
        for spin in range(self.nsppol):
            for ik, newk in enumerate(kfrac_coords):
                if dk1: der1 = dedk[spin, ik]
                if dk2: der2 = dedk2[spin, ik]
                new_eigens[spin, ik] = self.eval_sk(spin, newk, der1=der1, der2=der2)

        if self.verbose:
            print("Interpolation completed in %.3f (s)" % (time.time() - start))

        return dict2namedtuple(eigens=new_eigens, dedk=dedk, dedk2=dedk2)

    def interp_kpts_and_enforce_degs(self, kfrac_coords, ref_eigens, atol=1e-4):
        """
        Interpolate energies on an arbitrary set of k-points. Use `ref_eigens`
        to detect degeneracies and average the interpolated values in the degenerate subspace.
        """
        kfrac_coords = np.reshape(kfrac_coords, (-1, 3))
        new_nkpt = len(kfrac_coords)
        ref_eigens = np.reshape(ref_eigens, (self.nsppol, new_nkpt, self.nband))

        # Interpolate eigenvales.
        new_eigens = self.interp_kpts(kfrac_coords).eigens

        # Average interpolated values over degenerates bands.
        for spin in range(self.nsppol):
            for ik in range(new_nkpt):
                for dgbs in find_degs_sk(ref_eigens[spin, ik], atol):
                    if len(dgbs) == 1: continue
                    new_eigens[spin, ik, dgbs] = new_eigens[spin, ik, dgbs].sum() / len(dgbs)

        return dict2namedtuple(eigens=new_eigens, dedk=None, dedk2=None)

    def eval_sk(self, spin, kpt, der1=None, der2=None) -> np.ndarray:
        """
        Interpolate eigenvalues for all bands at a given (spin, k-point).
        Optionally compute gradients and Hessian matrices.

        Args:
            spin: Spin index.
            kpt: K-point in reduced coordinates.
            der1: If not None, ouput gradient is stored in der1[nband, 3].
            der2: If not None, output Hessian is der2[nband, 3, 3].

        Return:
            oeigs[nband]
        """
        # [NB, NR] x [NR]
        oeigs = np.matmul(self.coefs[spin], self.get_stark(kpt))
        if not self.iscomplexobj: oeigs = oeigs.real

        if der1 is not None:
            skr_dk1 = self.get_stark_dk1(kpt)
            for ii in range(3):
                value = np.matmul(self.coefs[spin, :, :], skr_dk1[ii])
                if not self.iscomplexobj: value = value.real
                der1[:, ii] = value

        if der2 is not None:
            skr_dk2 = self.get_stark_dk2(kpt)
            for jj in range(3):
                for ii in range(jj + 1):
                    value = np.matmul(self.coefs[spin, :, :], skr_dk2[ii,jj])
                    if not self.iscomplexobj: value = value.real
                    der2[:, ii, jj] = value
                    if ii != jj: der2[jj, ii] = der2[ii, jj]

        return oeigs

    #def eval_skb(self, spin, kpt, band, der1=None, der2=None):
    #    """
    #    Interpolate eigenvalues for a given (spin, k-point, band).

    #    Args:

    #    Return:
    #    """
    #    # Compute star function for this k-point (if not already in memory)
    #    if np.any(kpt != self.cached_kpt):
    #        self.cached_skr, self.cached_kpt = self.get_stark(kpt), kpt

    #    # Do not take complex conjugate.
    #    oeig = np.dot(self.coefs[spin, band], self.cached_skr)
    #    if not self.iscomplexobj: oeigs = oeigs.real
    #    if der1 is None and der2 is None: return oeig

    #    # TODO: Test Derivatives
    #    if der1 is not None:
    #        # Compute first-order derivatives.
    #        if np.any(kpt != self.cached_kpt_dk1):
    #            self.cached_skr_dk1, self.cached_kpt_dk1 = self.get_stark_dk1(kpt), kpt

    #        for ii in range(3):
    #            value = np.dot(self.coefs[spin, band], self.cached_skr_dk1[ii])
    #            if not self.iscomplexobj: value = value.real
    #            der1[ii] = value

    #        if der2 is None: return oeig, der1

    #    if der2 is not None:
    #        # Compute second-order derivatives.
    #        if np.any(kpt != self.cached_kpt_dk2):
    #            self.cached_skr_dk2, self.cached_kpt_dk2 = self.get_stark_dk2(kpt), kpt

    #        der2 = zero
    #        for jj in range(3):
    #            for ii in range(jj + 1):
    #                value = np.dot(self.coefs[spin, band], self.cached_skr_dk2[ii,jj])
    #                if not self.iscomplexobj: value = value.real
    #                der2[ii, jj] = value
    #                if ii != jj: der2[jj, ii] = der2[ii, jj]

    #    return oeig, der1, der2

    def get_stark(self, kpt) -> np.ndarray:
        """
        Return the star function for k-point `kpt`.

        Args:
            kpt: K-point in reduced coordinates.

        Return:
            complex array of shape [self.nr]
        """
        two_pi = 2.0 * np.pi
        skr = np.zeros(self.nr, dtype=complex)
        _np_exp = np.exp
        for omat in self.ptg_symrel:
            sk = two_pi * np.matmul(omat.T, kpt)
            skr += _np_exp(1.j * np.matmul(self.rpts, sk))
        skr /= self.ptg_nsym

        return skr

    def get_stark_dk1(self, kpt) -> np.ndarray:
        """
        Compute the 1st-order derivative of the star function wrt k

        Args:
            kpt: K-point in reduced coordinates.

        Return:
            complex array [3, self.nr]  with the derivative of the
            star function wrt k in reduced coordinates.
        """
        srk_dk1 = np.zeros((3, self.nr), dtype=complex)
        two_pi = 2.0 * np.pi
        rpts_t = self.rpts.T

        for omat in self.ptg_symrel:
            sk = two_pi * np.matmul(omat.T, kpt)
            exp_skr = np.exp(1.j * np.matmul(self.rpts, sk))
            for ir, rr in enumerate(self.rpts):
                srk_dk1[:, ir] += exp_skr[ir] * np.matmul(omat, rr)
            #for ir, or in enumerate(np.matmul(omat, rpts_t).T):
            #    srk_dk1[:, ir] += exp_skr[ir] * or

        srk_dk1 *= 1.j / self.ptg_nsym
        return srk_dk1

    def get_stark_dk2(self, kpt) -> np.ndarray:
        """
        Compute the 2nd-order derivatives of the star function wrt k.

        Args:
            kpt: K-point in reduced coordinates.

        Return:
            Complex numpy array of shape [3, 3, self.nr] with the 2nd-order derivatives
            of the star function wrt k in reduced coordinates.
        """
        srk_dk2 = np.zeros((3, 3, self.nr), dtype=complex)
        raise NotImplementedError()
        #work = zero
        #do isym=1,self.ptg_nsym
        #   sk = two_pi * matmul(transpose(self%ptg_symrel(:,:,isym)), kpt)
        #   do ir=1,self%nr
        #     sr = matmul(self%ptg_symrel(:,:,isym), self%rpts(:,ir))
        #     eiskr = exp(j_dpc * dot_product(sk, self%rpts(:,ir)))
        #     do jj=1,3
        #        do ii=1,jj
        #            work(ii,jj,ir) = work(ii,jj,ir) + eiskr * sr(ii) * sr(jj)

        #work = - work / self.ptg_nsym

        #do jj=1,3
        #   do ii=1,jj
        #       srk_dk2(:, ii, jj) = work(ii, jj, :)
        #       if (ii /= jj) srk_dk2(:,jj,ii) = work(:,ii,jj)

        return srk_dk2

    #def find_stationary_points(self, kmesh, bstart=None, bstop=None, is_shift=None)
    #    k = self.get_sampling(kmesh, is_shift)
    #    if bstart is None: bstart = self.nelect // 2 - 1
    #    if bstop is None: bstop = self.nelect // 2
    #    nb = bstop - bstart + 1
    #    results = []
    #    for ik_ibz, kpt in enumerate(k.ibz):
    #        vk_b = self.eval_dk1(kpt, bstart, bstop)
    #        bands = []
    #        for ib, v in enumerate(vk_b):
    #            if v < atol: bands.append(ib + bstart)
    #        if bands:
    #            #results.append()
    #            for band in bands:
    #                d2k_b = self.eval_dk2(kpt, band)

    #    return results

    def _find_rstar_gen(self, nrwant, rmax) -> tuple:
        """
        Find all lattice points generating the stars inside the supercell defined by `rmax`

        Args:
            nrwant: Number of star-functions required.
            rmax: numpy array with the maximum number of cells along the 3 reduced directions.

        Returns:
            tuple: (rpts, r2vals, ok)
        """
        msize = (2 * rmax + 1).prod()
        rtmp = np.empty((msize, 3), dtype=int)
        r2tmp = np.empty(msize)
        if self.verbose: print("rmax", rmax, "msize:", msize)

        start = time.time()
        for cnt, l in enumerate(itertools.product(range(-rmax[0], rmax[0] + 1),
                                                  range(-rmax[1], rmax[1] + 1),
                                                  range(-rmax[2], rmax[2] + 1))):
            rtmp[cnt] = l
            r2tmp[cnt] = np.dot(l, np.matmul(self.rmet, l))

        if self.verbose: print("gen points", time.time() - start)

        start = time.time()
        # Sort r2tmp and rtmp
        iperm = np.argsort(r2tmp)
        r2tmp = r2tmp[iperm]
        rtmp = rtmp[iperm]
        #return rtmp, r2tmp, True

        # Find shells
        r2sh = np.empty(msize, dtype=int)    # Correspondence between R and shell index.
        shlim = np.empty(msize, dtype=int)   # For each shell, the index of the initial G-vector.
        nsh = 1
        r2sh[0] = 0
        shlim[0] = 0
        r2_prev = 0.0
        for ir in range(1, msize):
            if abs(r2tmp[ir] - r2_prev) > r2tmp[ir] * 1e-8:
                r2_prev = r2tmp[ir]
                shlim[nsh] = ir
                nsh += 1
            r2sh[ir] = nsh
        shlim[nsh] = msize
        if self.verbose:
            print("nshells", nsh)
            print("shells", time.time() - start)

        # Find R-points generating the stars.
        # This is the most CPU critical part. I think it's difficult to do better than this without cython
        start = time.time()
        rgen = deque()

        for ish in range(nsh):
            ss, ee = shlim[ish], shlim[ish + 1]
            if ss + 1 == ee:
                rgen.append(rtmp[ss])
                continue

            if True:
                # This is faster
                rs = set([tuple(rtmp[ss])])
                for ir in range(ss + 1, ee):
                    rvec = rtmp[ir]
                    if all(tuple(np.matmul(rot, rvec)) not in rs for rot in self.ptg_symrel):
                        rs.add(tuple(rvec))
            else:
                rs = deque([rtmp[ss]])
                for ir in range(ss + 1, ee):
                    for rot in self.ptg_symrel:
                        rot_r = np.matmul(rot, rtmp[ir])
                        if any(np.all(rot_r == x) for x in rs): break
                    else:
                        # We have new R-point.
                        #print("Adding new point")
                        rs.append(rtmp[ir])

            #print(len(rs), rs)
            rgen.extend(rs)
        #print(rgen)
        if self.verbose: print("stars", time.time() - start)

        start = time.time()
        rgen = np.array(rgen, dtype=int)
        nstars = len(rgen)

        # Store rpts and compute ||R||**2.
        ok = nstars >= nrwant
        nr = min(nstars, nrwant)
        rpts = rgen[:nr].copy()
        r2vals = np.empty(nr)
        for ir in range(nr):
            r2vals[ir] = np.dot(rpts[ir], np.matmul(self.rmet, rpts[ir]))

        if self.verbose:
            print("r2max ", rpts[nr-1])
            print("end ", time.time() - start)
            if self.verbose > 10:
                print("nstars:", nstars)
                for r, r2 in zip(rpts, r2vals):
                    print(r, r2)

        return rpts, r2vals, ok


def extract_point_group(symrel, has_timrev) -> tuple:
    """
    Extract the point group rotations from the spacegroup. Add time-reversal
    if spatial inversion is not present and `has_timrev`.

    Return:
        (ptg_symrel, ptg_symrec) with rotations in real- and reciprocal-space.
    """
    nsym = len(symrel)
    tmp_nsym = 1
    work_symrel = np.empty((2*nsym, 3, 3), dtype=int)
    work_symrel[0] = symrel[0]

    for isym in range(1, nsym):
        found = any(np.all(symrel[isym] == work_symrel[search]) for search in range(tmp_nsym))
        if not found:
            work_symrel[tmp_nsym] = symrel[isym]
            tmp_nsym += 1

    inversion_3d = -np.eye(3, dtype=int)
    has_inversion = any(np.all(w == inversion_3d) for w in work_symrel[:tmp_nsym])

    # Now we know the symmetries of the point group.
    ptg_nsym = 2 * tmp_nsym if not has_inversion and has_timrev else tmp_nsym
    ptg_symrel = np.empty((ptg_nsym, 3, 3), dtype=int)
    ptg_symrec = np.empty((ptg_nsym, 3, 3), dtype=int)

    ptg_symrel[:tmp_nsym] = work_symrel[:tmp_nsym]
    for isym in range(tmp_nsym):
        ptg_symrec[isym] = np.linalg.inv(ptg_symrel[isym]).T

    if not has_inversion and has_timrev:
        # Add inversion.
        ptg_symrel[tmp_nsym:] = -work_symrel[:tmp_nsym]
        for isym in range(tmp_nsym, ptg_nsym):
            ptg_symrec[isym] = np.linalg.inv(ptg_symrel[isym]).T

    return ptg_symrel, ptg_symrec, has_inversion
