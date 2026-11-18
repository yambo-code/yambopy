# Copyright (c) 2025-2026, Ignacio M. Alliati, Myrta Grüning, Ashwin Pillai and Claudio Attaccalite
# All rights reserved.
#
# This file is part of the yambopy project
#
"""
Floquet analysis of real-time (yambo_nl / Lumen) simulations.

The real-time code propagates the occupied Bloch states |v_kn(t)> and, when
``FLOrder`` is set, stores their coefficients on the unperturbed Kohn-Sham basis

    C_in(t,k) = <u_ki | v_kn(t)>

at 2*N_F+1 equally spaced phases of one period (t_s = t_0 + s T/(2N_F+1)),
plus one record at t_0+T and one at t_0+4T (databases ndb.RT_V_bands and
ndb.RT_V_bands_K_section).

This module contains three parts.

Part I  - Direct monodromy formulation (recommended), class ``FloquetAnalysis``.
          All propagated states at a given k are treated as one subspace:

              G = C0^+ C0 ,   M = C0^+ CT ,   L = G^-1 M ,   L Gamma = Gamma diag(zeta_a)

              xi_ka = -arg(zeta_a)/T   (mod omega)

          The periodic part of the Floquet mode a and its Floquet-Kohn-Sham
          (Sambe) coefficients follow from

              d_kai(t)      = exp(+i xi_ka (t-t0)) sum_n C_in(t,k) Gamma_na
              d~_kai(eta)   = 1/(2N_F+1) sum_s exp(+i eta omega (t_s-t0)) d_kai(t_s)
              d_kai(t)      = sum_eta exp(-i eta omega (t-t0)) d~_kai(eta)

          i.e. psi_ka(t) = exp(-i xi_ka (t-t0)) sum_eta exp(-i eta omega (t-t0)) sum_i d~_kai(eta) u_ki
          (same sign convention as the paper and as the legacy VbPP class:
          the sideband eta carries the energy xi + eta*omega).
          The extraction is validated with: exact phase grid / T / 4T timing,
          condition numbers of G and Gamma, closure of the propagated subspace,
          semigroup consistency L_4 ~ L^4, off-diagonal weight of Gamma^-1 L_4 Gamma,
          Floquet-edge weight at |eta| = N_F and reconstruction of all stored records.

Part II - Berry-phase reconstruction of the polarization from the Floquet states
          and harmonic response chi^(n) (comparison with the real-time spectrum).

Part III- Legacy single-state API (VbPP, FLeigenvectors, findallk_qe, calc_rho,
          get_qebands_path, get_qebands_interpolate, plot_2D_kdist) kept unchanged
          for backward compatibility with the previous tutorial.

Notes
-----
* Only the *propagated occupied subspace* is represented: the number of Floquet
  modes per k equals the number of propagated valence bands. Conduction-band
  Floquet modes are not computed (they are not propagated by the RT code).
* Quasienergies are defined modulo hbar*omega. Branch shifts (gauges) only
  relabel the sidebands: d~(eta) -> d~(eta+m). 'compact', 'principal' and
  'ks' gauges are plotting/bookkeeping choices, not different physics.
* With respect to the stand-alone scripts FLQuasi.py / FLQ_HG.py the harmonic
  index has the opposite sign: Lambda_here[eta] = Lambda_scripts[-eta].

Theory: A. Pillai et al., "Floquet analysis of excitons in nonlinear response
functions" (in preparation).
"""
from yambopy import *
from yambopy.plot import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as pltcol
import sys
import os
import re
import glob
import json
import itertools
from datetime import datetime
from netCDF4 import Dataset
from yambopy.units import ha2ev, fs2aut
try:
    from yambopy.units import hbar_eVfs
except ImportError:  # older yambopy versions
    hbar_eVfs = 0.6582119569  # eV*fs

AUT2FS   = 1.0/fs2aut          # atomic unit of time -> fs
BOHR2ANG = 0.529177210903


class FloquetError(RuntimeError):
    """Raised when the Floquet extraction or one of its consistency checks fails."""
    pass


##################################################################################
### Small numerical helpers
##################################################################################

def _relnorm(delta, ref):
    """||delta|| / ||ref|| (Frobenius), robust to a zero reference."""
    den = float(np.linalg.norm(ref)); num = float(np.linalg.norm(delta))
    return num/den if den else (0.0 if not num else float('inf'))


def _numeric(var):
    return np.asarray(np.ma.asarray(var[:], dtype=float).filled(np.nan))


def _wrap(x, period):
    """Map x into (-period/2, period/2]."""
    return x - period*np.round(x/period)


##################################################################################
### PART I-a : reading the sampled Bloch states
##################################################################################

class FloquetSamples():
    """
    Container of the time-sampled propagated states written by yambo_nl/Lumen
    when FLOrder > 0.

    Attributes
    ----------
    times_au  : (nt,)                 time of each record [a.u.]
    records   : (nt,)                 original record index in ndb.RT_V_bands
    C         : (nt, nk, nbasis, nslot) complex, C[t,k,i,n] = <u_ki|v_kn(t)>
    band_lo, band_hi : KS basis window (NLBands, 1-based)
    nf        : Floquet order N_F (number of stored sidebands is 2*N_F+1)
    omega_au  : driving frequency [Ha]
    path      : directory the data were read from (or None)
    """
    def __init__(self, times_au, records, C, band_lo, band_hi, nf, omega_au, path=None):
        self.times_au = np.asarray(times_au, float)
        self.records  = np.asarray(records, int)
        self.C        = np.asarray(C, complex)
        self.band_lo  = int(band_lo)
        self.band_hi  = int(band_hi)
        self.nf       = int(abs(nf))
        self.omega_au = float(omega_au)
        self.path     = path
        if self.C.ndim != 4 or len(self.times_au) != self.C.shape[0]:
            raise FloquetError('Inconsistent time/record dimensions in the sampled states')
        if self.band_hi - self.band_lo + 1 != self.C.shape[2]:
            raise FloquetError('NLBands window and KS basis size are inconsistent')
        if self.nf < 1:
            raise FloquetError('FLOrder must be >= 1: no Floquet sampling found in the database')
        if not np.isfinite(self.omega_au) or self.omega_au <= 0:
            raise FloquetError('Invalid driving frequency')

    @property
    def period_au(self):
        return 2.0*np.pi/self.omega_au

    @property
    def nk(self):
        return self.C.shape[1]

    @property
    def basis_bands(self):
        return np.arange(self.band_lo, self.band_hi+1)

    @classmethod
    def from_db(cls, folder='.', calc='SAVE'):
        """Fast (vectorized) reader of ndb.RT_V_bands / ndb.RT_V_bands_K_section / ndb.Nonlinear."""
        path = os.path.join(folder, calc)
        main    = os.path.join(path, 'ndb.RT_V_bands')
        section = os.path.join(path, 'ndb.RT_V_bands_K_section')
        for f in (main, section):
            if not os.path.isfile(f):
                raise FloquetError(f'Missing database {f}')
        with Dataset(main) as db, Dataset(section) as sec:
            tv = db['IO_TIME_points']; times = _numeric(tv).ravel()
            unit = str(getattr(tv, 'units', '')).strip().lower()
            if unit in ('fs', 'femtosecond', 'femtoseconds'): times = times*fs2aut
            elif unit not in ('', 'au', 'a.u.', 'a.u', 'atomic units', 'atomic_units'):
                raise FloquetError(f'Unknown time unit in ndb.RT_V_bands: {unit}')
            var  = sec['V_bands']; dims = list(var.dimensions)
            wanted = ['V_bands_IO_time_step', 'n_sp_pol', 'QP_nk', 'nb_full', 'RT_nbands', 'complex']
            if len(dims) != 6 or set(dims) != set(wanted):
                raise FloquetError(f'Unsupported V_bands dimensions: {dims}')
            if var.shape[dims.index('n_sp_pol')] != 1:
                raise FloquetError('Spin-polarized V_bands are not supported yet')
            a = np.ma.asarray(var[:], dtype=float).filled(np.nan)
            a = a.take(0, axis=dims.index('n_sp_pol'))
            remaining = [d for d in dims if d != 'n_sp_pol']
            target = ['V_bands_IO_time_step', 'QP_nk', 'RT_nbands', 'nb_full', 'complex']
            a = a.transpose([remaining.index(d) for d in target])
            C = a[..., 0] + 1j*a[..., 1]
            bands = _numeric(db['RT_bands_kpts']).ravel()[:2].round().astype(int)
            nf = int(abs(round(_numeric(db['IO_Floquet_order']).ravel()[0])))
            nt = min(len(times), C.shape[0])
            if 'IO_TIME_steps_last_nsteps' in db.variables:
                nlast = int(_numeric(db['IO_TIME_steps_last_nsteps']).ravel()[2])
                if 0 < nlast < nt: nt = nlast
            times, C = times[:nt], C[:nt]
        omega = read_driving_frequency(path)
        records = np.arange(len(times))
        finite = np.isfinite(times) & np.isfinite(C).all(axis=(1, 2, 3))
        if not finite.all():
            print(f'[Floquet] non-finite records excluded: {records[~finite].tolist()}')
        times, C, records = times[finite], C[finite], records[finite]
        # duplicated records written after the 4T checkpoint are harmless: drop them
        last_required = 2*nf+2
        keep = np.ones(len(times), bool)
        for i in range(1, len(times)):
            if times[i] == times[i-1] and records[i] > last_required:
                if not np.allclose(C[i], C[i-1], rtol=1e-10, atol=1e-13):
                    raise FloquetError('Duplicate time stamp carries different wavefunctions')
                keep[i] = False
        times, C, records = times[keep], C[keep], records[keep]
        if len(times) < 3 or np.any(np.diff(times) <= 0):
            raise FloquetError('V_bands time stamps are not strictly increasing')
        return cls(times, records, C, bands[0], bands[1], nf, omega, path=path)

    @classmethod
    def from_vbandsdb(cls, vbdb):
        """Build the samples from an already loaded YamboVbandsDB object."""
        tv = np.array([np.asarray(v) for v in vbdb.tvecs])          # [t,k,nb_full,basis]
        C  = tv.transpose(0, 1, 3, 2)                                # [t,k,basis,nb_full]
        times = np.asarray(vbdb.times, float)
        records = np.arange(1, len(times)+1)                         # YamboVbandsDB skips record 0
        return cls(times, records, C, vbdb.basis_index[0], vbdb.basis_index[1],
                   vbdb.fl_order, vbdb.freq, path=getattr(vbdb, 'vb_path', None))

    def __str__(self):
        s  = "\n * * * Floquet samples (ndb.RT_V_bands) * * * \n\n"
        s += f"Path           : {self.path}\n"
        s += f"Records        : {len(self.times_au)} (original indices {self.records.tolist()})\n"
        s += f"Number of k    : {self.nk}\n"
        s += f"KS basis bands : {self.band_lo} - {self.band_hi}\n"
        s += f"Stored slots   : {self.C.shape[3]}\n"
        s += f"Floquet order  : {self.nf}  ({2*self.nf+1} phase samples)\n"
        s += f"Frequency      : {self.omega_au*ha2ev:.6f} eV\n"
        s += f"Period         : {self.period_au*AUT2FS:.6f} fs\n"
        return s


def read_driving_frequency(path):
    """Driving frequency [Ha] from ndb.Nonlinear."""
    with Dataset(os.path.join(path, 'ndb.Nonlinear')) as db:
        v = db['Field_Freq_1']; omega = float(_numeric(v).ravel()[0])
        unit = str(getattr(v, 'units', '')).strip().lower()
    if unit == 'ev': omega /= ha2ev
    elif unit not in ('', 'au', 'a.u.', 'a.u', 'ha', 'hartree', 'atomic units', 'atomic_units'):
        raise FloquetError(f'Unknown frequency unit: {unit}')
    return omega


def exact_samples(times, records, period, nf, tol=1e-6):
    """
    Locate the 2*N_F+1 phase samples of the first period and the T and 4T records,
    using the original record numbers (record 0 is the initial time and is not used).

    Returns i0, fit (array of indices), iT, i4, timing errors (in units of T).
    """
    ns  = 2*nf + 1
    loc = {int(r): i for i, r in enumerate(records)}
    missing = [r for r in range(1, ns+3) if r not in loc]
    if missing:
        raise FloquetError(f'Missing Floquet records {missing}: was the simulation long enough?')
    fit = np.array([loc[r] for r in range(1, ns+1)])
    i0, iT, i4 = loc[1], loc[ns+1], loc[ns+2]
    tau = (np.asarray(times) - times[i0])/period
    errors = {'phase_grid': float(np.max(abs(tau[fit] - np.arange(ns)/ns))),
              'T'         : float(abs(tau[iT] - 1)),
              'fourT'     : float(abs(tau[i4] - 4))}
    if max(errors.values()) > tol:
        raise FloquetError(f'Sampled times are not on the exact Floquet grid: {errors}')
    return i0, fit, iT, i4, errors


def active_slots(C0, rtol=1e-10):
    """Indices of the stored band slots that carry a propagated state."""
    norms = np.linalg.norm(C0, axis=1)              # [k, slot]
    mx = norms.max(axis=0)
    return np.flatnonzero(mx > rtol*max(float(mx.max()), 1.0))


##################################################################################
### PART I-b : direct monodromy formulation (pure numpy, k-vectorized)
##################################################################################

def diagonalize_monodromy(C0, CT, C4, period, cond_max=1e12):
    """
    One-period propagator in the subspace of the propagated states and its
    Floquet decomposition, for all k at once.

    C0, CT, C4 : (nk, nbasis, nstates) coefficients at t0, t0+T, t0+4T

    Returns a dict with
      Gram (G), L, L4, Gamma (metric-normalised eigenvectors, columns = modes),
      multipliers zeta, xi_au (principal quasienergy, real), imag_au (log|zeta|/T)
      and the diagnostic arrays (one value per k):
      gram_condition, eigenvector_condition, closure_T, closure_4T, semigroup_4T,
      predicted_4T, offdiag_4T, gram_unitarity, multiplier_modulus.
    """
    C0 = np.asarray(C0, complex); CT = np.asarray(CT, complex); C4 = np.asarray(C4, complex)
    C0h = np.conj(np.swapaxes(C0, 1, 2))
    G   = C0h @ C0
    gcond = np.linalg.cond(G)
    bad = ~np.isfinite(gcond) | (gcond > cond_max)
    if bad.any():
        raise FloquetError(f'Ill-conditioned Gram matrix at k = {(np.flatnonzero(bad)+1).tolist()[:10]}')
    L  = np.linalg.solve(G, C0h @ CT)
    L4 = np.linalg.solve(G, C0h @ C4)
    z, V = np.linalg.eig(L)
    if np.any(abs(z) < 1e-300):
        raise FloquetError('Zero monodromy eigenvalue')
    # normalise eigenvectors with the metric G:  Gamma^+ G Gamma = 1 on the diagonal
    norm = np.einsum('kia,kij,kja->ka', np.conj(V), G, V).real
    if np.any(norm <= 0):
        raise FloquetError('Non-positive eigenvector norm in the metric G')
    V = V/np.sqrt(norm)[:, None, :]
    vcond = np.linalg.cond(V)
    bad = ~np.isfinite(vcond) | (vcond > cond_max)
    if bad.any():
        raise FloquetError(f'Defective/ill-conditioned Floquet eigenvectors at k = {(np.flatnonzero(bad)+1).tolist()[:10]}')
    xi   = -np.angle(z)/period
    imag = np.log(abs(z))/period
    # diagnostics
    nk = len(C0)
    Lp4 = np.linalg.matrix_power(L, 4)
    rot = np.linalg.solve(V, L4 @ V)
    diag = np.einsum('kaa->ka', rot)
    chk = {}
    chk['gram_condition']        = gcond
    chk['eigenvector_condition'] = vcond
    chk['closure_T']     = np.array([_relnorm(C0[k] @ L[k] - CT[k], CT[k]) for k in range(nk)])
    chk['closure_4T']    = np.array([_relnorm(C0[k] @ L4[k] - C4[k], C4[k]) for k in range(nk)])
    chk['semigroup_4T']  = np.array([_relnorm(L4[k] - Lp4[k], L4[k]) for k in range(nk)])
    chk['predicted_4T']  = np.array([_relnorm(C0[k] @ Lp4[k] - C4[k], C4[k]) for k in range(nk)])
    chk['offdiag_4T']    = np.array([_relnorm(rot[k] - np.diag(diag[k]), diag[k]) for k in range(nk)])
    chk['gram_unitarity']= np.array([_relnorm(np.conj(L[k]).T @ G[k] @ L[k] - G[k], G[k]) for k in range(nk)])
    chk['multiplier_modulus'] = np.max(abs(abs(z) - 1.0), axis=1)
    out = {'Gram': G, 'L': L, 'L4': L4, 'Gamma': V, 'multipliers': z,
           'xi_au': xi, 'imag_au': imag}
    out.update(chk)
    return out


def assign_modes_to_bands(Gamma):
    """
    Associate each Floquet mode with one propagated (initial) band by maximising
    sum |Gamma_na|^2 (Hungarian algorithm). Returns perm[k, n] = mode of band n
    and the weight |Gamma_n,perm|^2 (1 = unambiguous, ~1/N near degeneracies).
    """
    from scipy.optimize import linear_sum_assignment
    nk, nst, _ = Gamma.shape
    perm = np.empty((nk, nst), int); weight = np.empty((nk, nst))
    for k in range(nk):
        w = abs(Gamma[k])**2
        w = w/np.maximum(w.sum(axis=0, keepdims=True), 1e-300)
        rows, cols = linear_sum_assignment(-w)
        perm[k, rows] = cols; weight[k, rows] = w[rows, cols]
    return perm, weight


def sambe_coefficients(Cfit, tau, Gamma, xi, omega, nf, shift='dominant'):
    """
    Floquet-Kohn-Sham (Sambe) coefficients by the finite DFT of the periodic part.

    Cfit  : (ns, nk, nbasis, nstates) coefficients at the phase samples
    tau   : (ns,) t_s - t0
    Gamma : (nk, nstates, nmodes);  xi : (nk, nmodes) quasienergies [Ha]
    shift : 'dominant' -> for each (k,a) the branch xi + m*omega is chosen so that the
                          sideband with the largest weight is eta = 0 (eta then counts
                          the photons exchanged with respect to the initial state);
            'compact'  -> branch with |m| <= nf minimising the edge weight at |eta| = nf,
                          then the second moment of the distribution, then |m|
                          (criterion of the stand-alone FLQuasi.py script; when the
                          edge weights are at noise level it is not unique);
            otherwise an integer array (nk, nmodes) of branch shifts m (0 = principal).

    Returns Lambda (nk, nmodes, 2nf+1, nbasis) with harmonics -nf..nf,
            the shift m used, the edge fraction and the spectral moment.
    """
    harm = np.arange(-nf, nf+1); ns = len(harm)
    tau = np.asarray(tau, float)
    if len(tau) != ns:
        raise FloquetError('Number of phase samples must be 2*N_F+1')
    modes = np.einsum('tkbi,kia->tkba', Cfit, Gamma, optimize=True)     # [t,k,basis,mode]
    Bplus = np.exp(+1j*np.outer(harm*omega, tau))                       # [eta, t]
    nk, nmode = xi.shape

    def _dft(m):
        e = xi + m*omega                                                 # (nk,nmode)
        D = modes*np.exp(1j*tau[:, None, None]*e[None, :, :])[:, :, None, :]
        coef = np.einsum('ht,tkba->kahb', Bplus, D, optimize=True)/ns     # [k,a,eta,basis]
        w = np.sum(abs(coef)**2, axis=3)
        tot = np.maximum(w.sum(axis=2), 1e-300)
        edge = (w[:, :, 0] + w[:, :, -1])/tot
        mom  = np.sum((harm/max(nf, 1))**2*w, axis=2)/tot
        return coef, edge, mom

    if isinstance(shift, str) and shift == 'dominant':
        c0, _, _ = _dft(np.zeros((nk, nmode)))
        shift = harm[np.argmax(np.sum(abs(c0)**2, axis=3), axis=2)]
    if not (isinstance(shift, str) or shift is None):
        shift = np.broadcast_to(np.asarray(shift, int), (nk, nmode)).copy()
        Lam = np.empty((nk, nmode, ns, Cfit.shape[2]), complex)
        edge = np.empty((nk, nmode)); mom = np.empty((nk, nmode))
        for m in np.unique(shift):
            c, e, mo = _dft(np.full((nk, nmode), m))
            sel = shift == m
            Lam[sel] = c[sel]; edge[sel] = e[sel]; mom[sel] = mo[sel]
        return Lam, shift, edge, mom

    cand = [_dft(np.full((nk, nmode), m)) for m in harm]
    E = np.stack([c[1] for c in cand]); M = np.stack([c[2] for c in cand])   # [m,k,a]
    best = np.empty((nk, nmode), int)
    absm = abs(harm)
    for k in range(nk):
        for a in range(nmode):
            best[k, a] = np.lexsort((absm, M[:, k, a], E[:, k, a]))[0]
    Lam  = np.empty((nk, nmode, ns, Cfit.shape[2]), complex)
    for i in range(ns):
        sel = best == i
        Lam[sel] = cand[i][0][sel]
    edge = np.take_along_axis(E, best[None], 0)[0]
    mom  = np.take_along_axis(M, best[None], 0)[0]
    return Lam, harm[best], edge, mom


def reconstruct_modes(Lambda, xi, omega, tau, harmonics):
    """Floquet modes psi_a(t0+tau) = exp(-i xi tau) sum_eta exp(-i eta omega tau) Lambda[eta]
       returned as [t, k, basis, mode]."""
    tau = np.atleast_1d(np.asarray(tau, float))
    B = np.exp(-1j*np.outer(tau, np.asarray(harmonics)*omega))               # [t,eta]
    D = np.einsum('th,kahb->tkba', B, Lambda, optimize=True)
    return D*np.exp(-1j*tau[:, None, None]*xi[None, :, :])[:, :, None, :]


def reconstruct_states(Lambda, xi, omega, tau, harmonics, Gamma):
    """Propagated states C(t0+tau) = modes(tau) Gamma^-1, returned as [t, k, basis, state]."""
    modes = reconstruct_modes(Lambda, xi, omega, tau, harmonics)
    return np.einsum('tkba,kaj->tkbj', modes, np.linalg.inv(Gamma), optimize=True)


def read_ks_energies(save='SAVE'):
    """KS eigenvalues [Ha] as (nk, nbands) from SAVE/ns.db1 (first spin channel)."""
    with Dataset(os.path.join(save, 'ns.db1')) as db:
        v = db['EIGENVALUES']; en = _numeric(v)
        unit = str(getattr(v, 'units', '')).strip().lower()
    if en.ndim == 3: en = en[0]
    if unit == 'ev': en = en/ha2ev
    return en


##################################################################################
### PART I-c : Brillouin-zone helpers (regular grid, periodic interpolation)
##################################################################################

def _periodic_tree(kred):
    from scipy.spatial import cKDTree
    k = np.mod(np.asarray(kred, float), 1.0)
    k[np.isclose(k, 1.0, atol=1e-12)] = 0.0
    return k, cKDTree(k, boxsize=1.0)


def _permutation_cycles(perm):
    n = len(perm); seen = np.zeros(n, bool); cycles = []
    for i in range(n):
        if seen[i]: continue
        cyc = []; j = i
        while not seen[j]:
            seen[j] = True; cyc.append(j); j = int(perm[j])
        if j != i: raise FloquetError('k-point neighbour map is not a permutation')
        cycles.append(cyc)
    return cycles


def _axis_neighbor(kred, axis, tol):
    k, tree = _periodic_tree(kred); nk = len(k)
    for N in sorted([d for d in range(2, nk+1) if nk % d == 0], reverse=True):
        delta = np.zeros(3); delta[axis] = 1.0/N
        dist, idx = tree.query(np.mod(k + delta, 1.0), k=1)
        idx = np.asarray(idx, int)
        if float(np.max(dist)) > tol or len(np.unique(idx)) != nk: continue
        cycles = _permutation_cycles(idx)
        if {len(c) for c in cycles} == {N}:
            return {'N': int(N), 'neighbor': idx, 'cycles': cycles}
    return {'N': 1, 'neighbor': np.arange(nk), 'cycles': [[i] for i in range(nk)]}


def infer_regular_grid(kred, tol=2e-5):
    """
    Detect the regular (Monkhorst-Pack-like) grid of a full-BZ k-point list.
    Returns shape (N1,N2,N3), per-axis neighbour information and the varying axes.
    """
    axes = {a: _axis_neighbor(kred, a, tol) for a in range(3)}
    shape = tuple(axes[a]['N'] for a in range(3))
    if int(np.prod(shape)) != len(kred):
        raise FloquetError(f'Could not detect a regular k grid (shape {shape}, nk={len(kred)})')
    return shape, axes, [a for a in range(3) if shape[a] > 1]


class FloquetBZ():
    """
    Full Brillouin-zone information needed to expand IBZ quantities, interpolate
    along paths and compute Berry phases. Built from the SAVE used by the
    real-time run (for a driven system: the FixSymm SAVE).
    """
    def __init__(self, save='SAVE', k_tol=2e-5):
        from yambopy.dbs.latticedb import YamboLatticeDB
        fname = os.path.join(save, 'ns.db1')
        lat  = YamboLatticeDB.from_db_file(filename=fname, Expand=False)
        full = YamboLatticeDB.from_db_file(filename=fname, Expand=True)
        self.lat_ibz = lat; self.lat_full = full
        self.kpoints_red      = np.asarray(lat.red_kpoints, float)
        self.kpoints_full_red = np.asarray(full.red_kpoints, float)
        self.bz_to_ibz        = np.asarray(full.kpoints_indexes, int)
        self.symmetry_indexes = np.asarray(full.symmetry_indexes, int)
        self.sym_car          = np.asarray(lat.sym_car, float)
        self.antiunitary      = np.asarray(lat.time_rev_list, bool)
        self.lattice_bohr     = np.asarray(lat.lat, float)
        self.rlat             = np.asarray(lat.rlat, float)        # without 2pi, Bohr^-1
        self.volume           = float(abs(lat.lat_vol))
        self.spin_factor      = float(getattr(lat, 'nspin', 2))
        self.nocc             = int(lat.nbandsv)
        self.shape, self.axes, self.varying = infer_regular_grid(self.kpoints_full_red, k_tol)
        self.ks_energies_au   = read_ks_energies(save)

    @property
    def nk_ibz(self): return len(self.kpoints_red)

    @property
    def nk_full(self): return len(self.kpoints_full_red)

    def expand(self, values):
        """IBZ -> full BZ (first axis)."""
        return np.asarray(values)[self.bz_to_ibz]

    def check_field_symmetries(self, direction):
        """The BZ expansion must use only unitary symmetries that leave the field invariant."""
        e = np.asarray(direction, float); e = e/np.linalg.norm(e)
        for isym in np.unique(self.symmetry_indexes):
            if self.antiunitary[isym]:
                raise FloquetError('BZ expansion uses time reversal: not allowed for a driven state (use a FixSymm SAVE)')
            if not np.allclose(self.sym_car[isym] @ e, e, atol=2e-6):
                raise FloquetError('BZ expansion uses a symmetry broken by the field: regenerate the SAVE with ypp FixSymm')

    def interpolate(self, values_ibz, points_red):
        """Periodic linear interpolation of an IBZ quantity on reduced-coordinate points."""
        from scipy.interpolate import RegularGridInterpolator
        k = self.kpoints_full_red; shape = np.asarray(self.shape); axes = self.varying
        values = self.expand(values_ibz)
        ref = k[0]; ix = np.rint((k - ref)*shape).astype(int) % shape
        if len(np.unique(ix, axis=0)) != len(k): raise FloquetError('Full-mesh index collision')
        cube = np.empty(tuple(shape) + values.shape[1:], dtype=values.dtype)
        cube[tuple(ix.T)] = values
        for a in sorted(set(range(3)) - set(axes), reverse=True): cube = np.take(cube, 0, axis=a)
        pads = [(0, 1)]*len(axes) + [(0, 0)]*(values.ndim-1)
        periodic = np.pad(cube, pads, mode='wrap')
        pts = np.asarray(points_red, float)
        for a in set(range(3)) - set(axes):
            if np.max(abs((pts[:, a] - ref[a] + .5) % 1 - .5)) > 1e-6:
                raise FloquetError('Path leaves the sampled k plane')
        query = np.mod(pts - ref, 1.0)[:, axes]*shape[axes]
        f = RegularGridInterpolator(tuple(np.arange(shape[a]+1) for a in axes), periodic, bounds_error=True)
        return f(query)

    def first_bz(self):
        """Wigner-Seitz cell (vertices, hull, plotting frame) in Cartesian 1/Angstrom."""
        from scipy.spatial import Voronoi, ConvexHull
        axes = np.asarray(self.varying, int); dim = len(axes)
        B = self.rlat[axes]*2*np.pi/BOHR2ANG
        if dim == 2:
            u = B[0]/np.linalg.norm(B[0]); n = np.cross(B[0], B[1]); n /= np.linalg.norm(n)
            frame = np.stack([u, np.cross(n, u)], axis=1)
        elif dim == 3: frame = np.eye(3)
        else: raise FloquetError('Only 2D/3D Brillouin zones are supported')
        b = B @ frame
        shifts = np.array(list(itertools.product(range(-2, 3), repeat=dim)))
        points = shifts @ b; origin = np.flatnonzero(np.all(shifts == 0, axis=1))[0]
        vor = Voronoi(points); reg = vor.regions[vor.point_region[origin]]
        vertices = vor.vertices[reg]; hull = ConvexHull(vertices)
        if dim == 2: vertices = vertices[hull.vertices]; hull = ConvexHull(vertices)
        return vertices, hull, frame

    def first_bz_values(self, values_ibz, span=2):
        """Full-BZ points folded in the first BZ (Cartesian, plotting frame) and their values."""
        axes = np.asarray(self.varying, int)
        bz, hull, frame = self.first_bz(); dim = len(axes)
        vals = self.expand(values_ibz)
        shifts = np.zeros(((2*span+1)**dim, 3))
        shifts[:, axes] = np.array(list(itertools.product(range(-span, span+1), repeat=dim)))
        k = self.kpoints_full_red - np.floor(self.kpoints_full_red + .5)
        pts = (k[None] + shifts[:, None]).reshape(-1, 3) @ self.rlat*2*np.pi/BOHR2ANG @ frame
        allv = np.tile(vals, len(shifts))
        inside = np.all(pts @ hull.equations[:, :dim].T + hull.equations[:, dim] <= 1e-8, axis=1)
        return pts, allv, inside, bz, hull, frame


##################################################################################
### PART I-d : the FloquetAnalysis class
##################################################################################

class FloquetAnalysis():
    """
    Direct-monodromy Floquet analysis of the propagated occupied subspace,
    for all k-points of a single-frequency yambo_nl run with FLOrder > 0.

    Usage
    -----
    >>> fa = FloquetAnalysis(folder='.', calc='FLQ', save='SAVE')    # reads ndb.RT_V_bands*
    >>> print(fa)                       # summary + consistency checks
    >>> xi = fa.quasienergies(gauge='ks')        # (nk, nmodes) in eV
    >>> eta, d = fa.sambe(kpt=6, band=4)         # F-KS coefficients d~_kai(eta)
    >>> fa.plot_sambe(6, 4); fa.plot_periodic_part(6, 4, basis_band=5)

    An already loaded YamboVbandsDB can be passed instead: FloquetAnalysis(vbdb=mydb).

    Parameters
    ----------
    vbdb / samples : YamboVbandsDB or FloquetSamples (alternative to folder/calc)
    folder, calc   : location of ndb.RT_V_bands, ndb.RT_V_bands_K_section, ndb.Nonlinear
    save           : SAVE folder of the run (FixSymm SAVE for a driven system); used for
                     KS energies, the number of occupied bands and the BZ. Can be None.
    gauge          : branch used for the Sambe coefficients:
                     'dominant' (default) the largest sideband is eta=0, so eta counts the
                                photons exchanged with the initial band (eta=+-2 <-> SHG);
                     'compact'  minimal weight at |eta|=N_F (as in FLQuasi.py);
                     'principal' xi in (-omega/2, omega/2]
    time_tol       : tolerance on the Floquet time grid (units of T)
    tol            : tolerance on closure/semigroup/reconstruction tests
    strict         : raise FloquetError if a test fails (default), otherwise only warn
    """

    def __init__(self, vbdb=None, folder='.', calc=None, save='SAVE', samples=None,
                 gauge='dominant', time_tol=1e-6, tol=1e-6, cond_max=1e12,
                 active_rtol=1e-10, strict=True, verbose=True):
        if samples is not None: self.samples = samples
        elif vbdb is not None:  self.samples = FloquetSamples.from_vbandsdb(vbdb)
        elif calc is not None:  self.samples = FloquetSamples.from_db(folder, calc)
        else: raise ValueError('Provide vbdb, samples or calc')
        self.folder, self.calc = folder, calc
        self.save = save if (save is not None and os.path.isfile(os.path.join(save, 'ns.db1'))) else None
        if save is not None and self.save is None and verbose:
            print(f'[Floquet] {save}/ns.db1 not found: KS energies not available')
        self.ks_au = read_ks_energies(self.save) if self.save else None
        self.nocc = None
        if self.save:
            from yambopy.dbs.latticedb import YamboLatticeDB
            self.nocc = int(YamboLatticeDB.from_db_file(filename=os.path.join(self.save, 'ns.db1'),
                                                        Expand=False).nbandsv)
            if self.ks_au.shape[0] != self.samples.nk:
                raise FloquetError(f'SAVE has {self.ks_au.shape[0]} k-points, V_bands has {self.samples.nk}: use the SAVE of the RT run')
        self.gauge, self.time_tol, self.tol, self.cond_max = gauge, time_tol, tol, cond_max
        self.active_rtol, self.strict, self.verbose = active_rtol, strict, verbose
        self.bz = None
        self.run()

    # ------------------------------------------------------------------ core
    def run(self):
        s = self.samples
        self.omega_au = s.omega_au; self.period_au = s.period_au; self.nf = s.nf
        i0, fit, iT, i4, self.timing_errors = exact_samples(s.times_au, s.records, s.period_au, s.nf, self.time_tol)
        self.i0, self.fit, self.iT, self.i4 = i0, fit, iT, i4
        self.t0_au = float(s.times_au[i0])
        slots = active_slots(s.C[i0], self.active_rtol)
        if len(slots) == 0: raise FloquetError('No propagated state found in V_bands')
        if slots[0] == s.band_lo-1:   bands = slots + 1             # absolute band slots
        elif slots[0] == 0:           bands = s.band_lo + slots      # local slots
        else: raise FloquetError(f'Cannot map stored slots {slots.tolist()} to KS bands')
        if np.any(np.diff(bands) != 1) or bands[0] < s.band_lo or bands[-1] > s.band_hi:
            raise FloquetError(f'Propagated bands {bands.tolist()} not contiguous or outside NLBands')
        if self.nocc is not None and bands[-1] != self.nocc:
            raise FloquetError(f'Propagated bands {bands.tolist()} do not end at the last occupied band {self.nocc}')
        self.active_slots, self.active_bands = slots, bands
        self.basis_bands = s.basis_bands
        C = s.C[..., slots]
        self.C = C
        d = diagonalize_monodromy(C[i0], C[iT], C[i4], self.period_au, self.cond_max)
        # order the modes following the propagated bands they originate from
        perm, weight = assign_modes_to_bands(d['Gamma'])
        take = lambda a: np.take_along_axis(a, perm, axis=1)
        self.Gamma       = np.take_along_axis(d['Gamma'], perm[:, None, :], axis=2)
        self.multipliers = take(d['multipliers'])
        self.xi_principal_au = take(d['xi_au'])
        self.imag_au     = take(d['imag_au'])
        self.mode_weight = weight
        self.Gram, self.L, self.L4 = d['Gram'], d['L'], d['L4']
        self.diagnostics = {k: d[k] for k in ('gram_condition', 'eigenvector_condition', 'closure_T',
                            'closure_4T', 'semigroup_4T', 'predicted_4T', 'offdiag_4T',
                            'gram_unitarity', 'multiplier_modulus')}
        if self.gauge not in ('dominant', 'compact', 'principal'):
            raise ValueError("gauge must be 'dominant', 'compact' or 'principal'")
        shift = 0 if self.gauge == 'principal' else self.gauge
        tau = s.times_au[fit] - self.t0_au
        self.Lambda, self.gauge_shift, edge, mom = sambe_coefficients(
            C[fit], tau, self.Gamma, self.xi_principal_au, self.omega_au, self.nf, shift=shift)
        self.harmonics = np.arange(-self.nf, self.nf+1)
        self.xi_au = self.xi_principal_au + self.gauge_shift*self.omega_au
        self.diagnostics['edge_fraction'] = edge
        self.diagnostics['spectral_moment'] = mom
        # reconstruction of every stored record after t0 (phase samples, T and 4T)
        idx = np.flatnonzero(s.times_au >= self.t0_au)
        rec = self.reconstruct_states(s.times_au[idx] - self.t0_au)
        self.reconstruction_records = s.records[idx]
        self.reconstruction_errors = np.array([_relnorm(a-b, b) for a, b in zip(rec, C[idx])])
        self.passed = self._check()

    def _check(self):
        failed = []
        for key in ('closure_T', 'closure_4T', 'semigroup_4T', 'predicted_4T'):
            if self.diagnostics[key].max() > self.tol: failed.append(key)
        if self.reconstruction_errors.max() > self.tol: failed.append('reconstruction')
        if failed:
            msg = f'Floquet consistency tests failed (tol={self.tol:g}): {failed}'
            if self.strict: raise FloquetError(msg)
            print('[WARNING] '+msg)
        if self.verbose and self.diagnostics['edge_fraction'].max() > 1e-6:
            print(f"[WARNING] weight at |eta|=N_F up to {self.diagnostics['edge_fraction'].max():.2e}: "
                  "increase FLOrder")
        return not failed

    # ------------------------------------------------------------- utilities
    def _k(self, kpt):
        if not 1 <= kpt <= self.samples.nk: raise ValueError(f'kpt must be in 1..{self.samples.nk}')
        return kpt-1

    def _mode(self, band):
        """Index of the Floquet mode associated with the propagated KS band (1-based)."""
        if band is None: return len(self.active_bands)-1
        bl = list(self.active_bands)
        if band not in bl: raise ValueError(f'band must be one of the propagated bands {bl}')
        return bl.index(band)

    def _basis(self, basis_band):
        bl = list(self.basis_bands)
        if basis_band not in bl: raise ValueError(f'basis_band must be one of {bl}')
        return bl.index(basis_band)

    @property
    def nmodes(self): return len(self.active_bands)

    @property
    def nk(self): return self.samples.nk

    @property
    def ks_offset_au(self):
        """
        Difference between the energy reference of the RT run and the KS energies in
        SAVE, modulo omega. yambo_nl puts the valence-band maximum at zero, so the
        offset is -E_VBM (mod omega). This is a reference shift only: the k-dependent
        differences xi - E_KS (Stark/Bloch-Siegert shifts) are left untouched.
        See median_offset_au for the value estimated from the data.
        """
        if self.ks_au is None: return 0.0
        vbm = self.ks_au[:, self.active_bands[-1]-1].max()
        return float(_wrap(-vbm, self.omega_au))

    @property
    def median_offset_au(self):
        """Median over k and modes of (xi - E_KS) mod omega: should be close to ks_offset_au,
           the difference being the average field-induced shift of the quasienergies."""
        if self.ks_au is None: return 0.0
        res = _wrap(self.xi_principal_au - self.ks_au[:, self.active_bands-1], self.omega_au)
        return float(np.median(res))

    def quasienergies(self, gauge='ks', align=True, units='eV'):
        """
        Quasienergies (nk, nmodes); mode a is associated with band active_bands[a].
          gauge = 'principal' : -arg(zeta)/T in (-omega/2, omega/2]
                  'sambe'     : branch used for the Sambe coefficients Lambda (gauge of the class)
                  'ks'        : branch closest to the KS energy of the associated band;
                                with align=True the global offset ks_offset_au is removed so
                                that xi can be directly compared with the KS energies of SAVE.
        """
        if gauge == 'principal': xi = self.xi_principal_au.copy()
        elif gauge in ('sambe', self.gauge): xi = self.xi_au.copy()
        elif gauge == 'ks':
            if self.ks_au is None: raise FloquetError("gauge='ks' needs the KS energies (save=...)")
            off = self.ks_offset_au if align else 0.0
            eks = self.ks_au[:, self.active_bands-1]
            n = np.rint((eks + off - self.xi_principal_au)/self.omega_au)
            xi = self.xi_principal_au + n*self.omega_au - off
        else: raise ValueError("gauge must be 'principal', 'sambe' or 'ks'")
        return xi*ha2ev if units == 'eV' else xi

    def ks_energies(self, bands=None, units='eV'):
        """KS energies (nk, nbands) of SAVE for the given 1-based bands (default: propagated bands)."""
        if self.ks_au is None: raise FloquetError('KS energies not available (save=None)')
        bands = self.active_bands if bands is None else np.atleast_1d(bands)
        e = self.ks_au[:, np.asarray(bands)-1]
        return e*ha2ev if units == 'eV' else e

    def sambe(self, kpt, band=None):
        """Harmonics eta and F-KS coefficients d~_kai(eta), array (2N_F+1, nbasis), for the
           Floquet mode associated with the propagated band `band` at k-point `kpt` (1-based)."""
        return self.harmonics, self.Lambda[self._k(kpt), self._mode(band)]

    def sideband_weight(self, eta, basis_bands=None, bands=None, quantity='weight'):
        """
        Weight of sideband eta for every k:  sum_{a in bands} sum_{i in basis_bands} |d~_kai(eta)|^2
        (quantity='amplitude' returns the square root).
        Default: complete sums over modes and basis bands, which are invariant under
        rotations within degenerate subspaces (mode/band-resolved maps may not be).
        """
        ih = np.flatnonzero(self.harmonics == eta)
        if len(ih) != 1: raise ValueError(f'eta must be in {self.harmonics.tolist()}')
        ib = slice(None) if basis_bands is None else [self._basis(b) for b in np.atleast_1d(basis_bands)]
        im = slice(None) if bands is None else [self._mode(b) for b in np.atleast_1d(bands)]
        c = self.Lambda[:, :, ih[0], :][:, im][:, :, ib]
        w = np.sum(abs(c)**2, axis=(1, 2))
        return np.sqrt(w) if quantity == 'amplitude' else w

    def reconstruct_modes(self, tau):
        """Floquet modes at t0+tau, array [t, k, basis, mode]."""
        return reconstruct_modes(self.Lambda, self.xi_au, self.omega_au, tau, self.harmonics)

    def reconstruct_states(self, tau):
        """Propagated states at t0+tau rebuilt from the Floquet modes, [t, k, basis, state]."""
        return reconstruct_states(self.Lambda, self.xi_au, self.omega_au, tau, self.harmonics, self.Gamma)

    def periodic_part(self, kpt, band, tau):
        """d_kai(t0+tau) = sum_eta exp(-i eta omega tau) d~_kai(eta), array [t, basis]."""
        k, a = self._k(kpt), self._mode(band)
        B = np.exp(-1j*np.outer(np.atleast_1d(tau), self.harmonics*self.omega_au))
        return B @ self.Lambda[k, a]

    def raw_periodic_part(self, kpt, band):
        """Periodic part obtained directly from the stored records:
           exp(+i xi (t-t0)) sum_n C_in(t) Gamma_na. Returns (tau, [t, basis])."""
        k, a = self._k(kpt), self._mode(band)
        idx = np.flatnonzero(self.samples.times_au >= self.t0_au)
        tau = self.samples.times_au[idx] - self.t0_au
        d = (self.C[idx, k] @ self.Gamma[k][:, a])*np.exp(1j*self.xi_au[k, a]*tau)[:, None]
        return tau, d

    # --------------------------------------------------------------- reports
    def check_summary(self):
        """min/mean/max of every diagnostic."""
        out = {}
        for key, v in self.diagnostics.items():
            out[key] = (float(np.min(v)), float(np.mean(v)), float(np.max(v)))
        r = self.reconstruction_errors
        out['reconstruction'] = (float(r.min()), float(r.mean()), float(r.max()))
        return out

    def __str__(self):
        s  = "\n * * * Floquet analysis (direct monodromy) * * * \n"
        s += str(self.samples)
        s += f"Propagated bands  : {self.active_bands.tolist()}\n"
        s += f"t0                : {self.t0_au*AUT2FS:.6f} fs\n"
        s += f"Sambe gauge       : {self.gauge}\n"
        if self.ks_au is not None:
            s += f"KS/RT offset      : {self.ks_offset_au*ha2ev:+.6f} eV (-VBM mod omega); from data: {self.median_offset_au*ha2ev:+.6f} eV\n"
        s += f"Timing errors / T : " + ", ".join(f"{k}={v:.1e}" for k, v in self.timing_errors.items()) + "\n"
        s += "\n Consistency checks (min / mean / max over k)\n"
        for key, (a, b, c) in self.check_summary().items():
            s += f"   {key:22s} {a:.3e}  {b:.3e}  {c:.3e}\n"
        s += f"\n All tests passed  : {self.passed}\n"
        return s

    def report(self, filename=None):
        """Per-k table of KS energies, quasienergies (KS gauge), assignment weight and checks."""
        lines = ['#  k  band   E_KS[eV]      xi[eV]      xi-E_KS[eV]  weight   edge      closure_T']
        xi = self.quasienergies('ks') if self.ks_au is not None else self.quasienergies('sambe')
        eks = self.ks_energies() if self.ks_au is not None else np.full_like(xi, np.nan)
        for k in range(self.nk):
            for a, b in enumerate(self.active_bands):
                lines.append(f'{k+1:4d} {b:4d} {eks[k,a]:12.6f} {xi[k,a]:12.6f} {xi[k,a]-eks[k,a]:12.4e} '
                             f'{self.mode_weight[k,a]:7.4f} {self.diagnostics["edge_fraction"][k,a]:9.2e} '
                             f'{self.diagnostics["closure_T"][k]:9.2e}')
        text = '\n'.join(lines) + '\n'
        if filename:
            with open(filename, 'w') as f: f.write(text)
        else: print(text)

    def to_npz(self, filename):
        """Save quasienergies, Sambe coefficients and diagnostics in a numpy archive."""
        data = dict(omega_au=self.omega_au, period_au=self.period_au, t0_au=self.t0_au, nf=self.nf,
                    harmonics=self.harmonics, active_bands=self.active_bands, basis_bands=self.basis_bands,
                    xi_principal_au=self.xi_principal_au, xi_au=self.xi_au, gauge_shift=self.gauge_shift,
                    imag_au=self.imag_au, multipliers=self.multipliers, Gamma=self.Gamma, Lambda=self.Lambda,
                    mode_weight=self.mode_weight, times_au=self.samples.times_au, records=self.samples.records,
                    C=self.C, reconstruction_errors=self.reconstruction_errors,
                    convention='psi=exp(-i xi tau) sum_eta exp(-i eta omega tau) Lambda[eta]')
        for k, v in self.diagnostics.items(): data['check_'+k] = v
        if self.ks_au is not None: data['ks_energies_au'] = self.ks_au
        np.savez_compressed(filename, **data)

    # ----------------------------------------------------------------- plots
    def plot_sambe(self, kpt, band=None, ax=None, vmin=None, cmap='YlOrBr', annotate=True, savefig=None):
        """Map of |d~_kai(eta)|: KS basis band (rows) vs Floquet harmonic eta (columns).
           vmin: lower limit of the log colour scale (default 1e-12 of the largest value)."""
        eta, d = self.sambe(kpt, band)
        if vmin is None: vmin = max(abs(d).max()*1e-12, 1e-300)
        amp = np.maximum(abs(d.T), vmin)
        if ax is None: fig, ax = plt.subplots(figsize=(8, 3.5))
        else: fig = ax.figure
        im = ax.imshow(amp, origin='lower', aspect='auto', cmap=cmap,
                       norm=pltcol.LogNorm(vmin=vmin, vmax=max(1.0, amp.max())))
        ax.set_xticks(range(len(eta))); ax.set_xticklabels(eta)
        ax.set_yticks(range(len(self.basis_bands))); ax.set_yticklabels(self.basis_bands)
        for i in range(amp.shape[0]):
            for j in range(amp.shape[1]):
                if annotate and amp[i, j] > 1e3*vmin:
                    ax.text(j, i, f'{amp[i,j]:.1e}', ha='center', va='center', fontsize=7)
        ax.set_xlabel(r'Floquet harmonic $\eta$'); ax.set_ylabel('KS basis band $i$')
        b = self.active_bands[self._mode(band)]
        ax.set_title(rf'$|\tilde d_{{kai}}(\eta)|$   k={kpt}, mode of band {b}, '
                     rf'$\xi$={self.xi_au[kpt-1,self._mode(band)]*ha2ev:.4f} eV ({self.gauge} gauge)')
        fig.colorbar(im, ax=ax)
        if savefig: fig.savefig(savefig, bbox_inches='tight', dpi=200)
        return ax

    def plot_periodic_part(self, kpt, band, basis_band, nperiods=4, npts=1201, savefig=None):
        """Re/Im of the periodic part d_kai(t): Fourier reconstruction (line) vs stored records (dots)."""
        ib = self._basis(basis_band)
        tau_raw, draw = self.raw_periodic_part(kpt, band)
        tau = np.linspace(0, nperiods*self.period_au, npts)
        d = self.periodic_part(kpt, band, tau)[:, ib]
        fig, axes = plt.subplots(2, sharex=True, figsize=(8, 5.5))
        for ax, f, lab in zip(axes, (np.real, np.imag), ('Re', 'Im')):
            ax.plot((tau+self.t0_au)*AUT2FS, f(d), color='tab:blue', label='Floquet reconstruction')
            ax.plot((tau_raw+self.t0_au)*AUT2FS, f(draw[:, ib]), 'o', color='tab:olive', ms=8,
                    label='RT records')
            for p in (1, 4): ax.axvline((self.t0_au+p*self.period_au)*AUT2FS, ls='--', lw=.7, color='0.5')
            ax.set_ylabel(rf'{lab} $d_{{k a i}}(t)$')
        axes[1].set_xlabel('Time (fs)'); axes[0].legend(fontsize=8)
        axes[0].set_title(f'Periodic part: k={kpt}, mode of band {band}, KS basis band {basis_band}')
        if savefig: fig.savefig(savefig, bbox_inches='tight', dpi=200)
        return axes

    def plot_propagated_state(self, kpt, band, basis_band, nperiods=4, npts=1201, savefig=None):
        """Re/Im/|.| of the propagated coefficient C_in(t) rebuilt from all Floquet modes vs RT records."""
        k, j, ib = self._k(kpt), self._mode(band), self._basis(basis_band)
        tau = np.linspace(0, nperiods*self.period_au, npts)
        c = self.reconstruct_states(tau)[:, k, ib, j]
        idx = np.flatnonzero(self.samples.times_au >= self.t0_au)
        traw = self.samples.times_au[idx]; craw = self.C[idx, k, ib, j]
        fig, axes = plt.subplots(3, sharex=True, figsize=(8, 7))
        for ax, f, lab in zip(axes, (np.real, np.imag, np.abs), ('Re', 'Im', '|.|')):
            ax.plot((tau+self.t0_au)*AUT2FS, f(c), color='tab:blue', label='Floquet reconstruction')
            ax.plot(traw*AUT2FS, f(craw), 'o', color='tab:olive', ms=8, label='RT records')
            ax.set_ylabel(rf'{lab} $C_{{in}}(t)$')
        axes[2].set_xlabel('Time (fs)'); axes[0].legend(fontsize=8)
        axes[0].set_title(f'Propagated state: k={kpt}, band n={band}, KS basis band i={basis_band}')
        if savefig: fig.savefig(savefig, bbox_inches='tight', dpi=200)
        return axes

    def plot_checks(self, savefig=None):
        """Diagnostics as a function of the k index (log scale)."""
        fig, ax = plt.subplots(figsize=(8, 4.5))
        k = np.arange(1, self.nk+1)
        for key in ('closure_T', 'closure_4T', 'semigroup_4T', 'offdiag_4T', 'gram_unitarity'):
            ax.semilogy(k, np.maximum(self.diagnostics[key], 1e-20), '.', ms=3, label=key)
        ax.semilogy(k, np.maximum(self.diagnostics['edge_fraction'].max(axis=1), 1e-30), '.', ms=3,
                    label=r'edge weight $|\eta|=N_F$')
        ax.axhline(self.tol, color='k', ls='--', lw=.8)
        ax.set_xlabel('k index'); ax.set_ylabel('relative error'); ax.legend(fontsize=7, ncol=2)
        if savefig: fig.savefig(savefig, bbox_inches='tight', dpi=200)
        return ax

    def _bz(self, bz):
        if bz is not None: self.bz = bz
        if self.bz is None:
            if self.save is None: raise FloquetError('A FloquetBZ (or save=...) is needed')
            self.bz = FloquetBZ(self.save)
        if self.bz.nk_ibz != self.nk:
            raise FloquetError('FloquetBZ and V_bands have a different number of k-points')
        return self.bz

    def bands_on_path(self, path, npoints=30, bz=None, ks_bands=None, align=True):
        """
        Floquet quasienergies (KS-connected branch) and KS energies along a path.

        The Floquet multipliers exp(-i xi T), which do not depend on the branch, are
        interpolated on the periodic full-BZ grid; the branch closest to the
        interpolated KS energy of the associated band is then selected.

        path : qepy Path, or list of [[k1,k2,k3],'label'] in reduced coordinates
        Returns dict with distance (1/Ang), ticks, labels, xi (eV), ks (eV), points.
        """
        bz = self._bz(bz)
        if hasattr(path, 'kpoints'):
            verts, labels = np.asarray(path.kpoints, float), list(path.klabels)
            nseg = list(path.intervals) if npoints is None else [npoints]*(len(verts)-1)
        else:
            verts = np.array([p[0] for p in path], float); labels = [p[1] for p in path]
            nseg = [npoints or 30]*(len(verts)-1)
        B = bz.rlat*2*np.pi/BOHR2ANG
        pts, dist, ticks, off = [], [], [0.0], 0.0
        for (a, b), n in zip(zip(verts[:-1], verts[1:]), nseg):
            seg = np.linspace(a, b, n, endpoint=False); length = float(np.linalg.norm((b-a) @ B))
            pts.append(seg); dist.extend(off + np.arange(n)*length/n); off += length; ticks.append(off)
        pts = np.vstack(pts + [verts[-1:]]); dist = np.r_[dist, off]
        ks_bands = np.arange(1, self.active_bands[-1]+3) if ks_bands is None else np.asarray(ks_bands)
        ks = bz.interpolate(self.ks_au[:, ks_bands-1], pts)
        ks_assoc = bz.interpolate(self.ks_au[:, self.active_bands-1], pts)
        z = bz.interpolate(np.exp(-1j*self.xi_principal_au*self.period_au), pts)
        if np.min(abs(z)) < 0.1:
            raise FloquetError('Interpolated Floquet multipliers vanish: refine the k grid')
        principal = -np.angle(z)/self.period_au
        off_ks = self.ks_offset_au if align else 0.0
        n = np.rint((ks_assoc + off_ks - principal)/self.omega_au)
        xi = principal + n*self.omega_au - off_ks
        return {'distance': dist, 'ticks': np.array(ticks), 'labels': labels, 'points': pts,
                'xi': xi*ha2ev, 'ks': ks*ha2ev, 'ks_bands': ks_bands, 'ks_assoc': ks_assoc*ha2ev}

    def plot_bands(self, path, npoints=30, bz=None, ks_bands=None, replicas=1, align=True,
                   zero='vbm', ax=None, ylim=None, savefig=None):
        """Floquet quasienergies of the propagated modes (and +-m*hbar*omega replicas) vs KS bands."""
        r = self.bands_on_path(path, npoints, bz, ks_bands, align)
        e0 = self.ks_au[:, self.active_bands[-1]-1].max()*ha2ev if zero == 'vbm' else float(zero)
        w = self.omega_au*ha2ev
        if ax is None: fig, ax = plt.subplots(figsize=(6, 6))
        else: fig = ax.figure
        for j in range(r['ks'].shape[1]):
            ax.plot(r['distance'], r['ks'][:, j]-e0, color='tab:red', lw=1.2, label='KS' if j == 0 else None)
        for a in range(self.nmodes):
            for m in range(-replicas, replicas+1):
                ax.plot(r['distance'], r['xi'][:, a]+m*w-e0, '--', color='tab:blue', lw=1.,
                        alpha=1.0 if m == 0 else .5,
                        label=(r'Floquet $\xi$ (+ $m\hbar\omega$)' if a == 0 and m == 0 else None))
        for t in r['ticks']: ax.axvline(t, color='0.5', lw=.6, ls=':')
        ax.set_xticks(r['ticks']); ax.set_xticklabels(r['labels']); ax.set_xlim(r['distance'][0], r['distance'][-1])
        ax.set_ylabel('Energy (eV)' + (' w.r.t. VBM' if zero == 'vbm' else ''))
        ax.set_title(rf'$\hbar\omega$ = {w:.3f} eV')
        if ylim: ax.set_ylim(ylim)
        ax.legend(fontsize=8, loc='lower left')
        if savefig: fig.savefig(savefig, bbox_inches='tight', dpi=200)
        return ax

    def plot_sideband_bz(self, eta, basis_bands=None, bands=None, quantity='weight', bz=None,
                         grid=250, path=None, ax=None, cmap='viridis', show_kpoints=False, savefig=None):
        """
        Map of the sideband weight sum |d~_kai(eta)|^2 in the first Brillouin zone (2D systems).
        It is a Floquet-coefficient weight, not an exciton wavefunction.
        """
        from scipy.interpolate import LinearNDInterpolator
        bz = self._bz(bz)
        vals = self.sideband_weight(eta, basis_bands, bands, quantity)
        pts, allv, inside, cell, hull, frame = bz.first_bz_values(vals)
        if len(bz.varying) != 2: raise FloquetError('plot_sideband_bz is implemented for 2D systems')
        if ax is None: fig, ax = plt.subplots(figsize=(6.5, 5.5))
        else: fig = ax.figure
        xmin, ymin = cell.min(axis=0); xmax, ymax = cell.max(axis=0)
        X, Y = np.meshgrid(np.linspace(xmin, xmax, grid), np.linspace(ymin, ymax, grid))
        tgt = np.c_[X.ravel(), Y.ravel()]
        mask = np.all(tgt @ hull.equations[:, :2].T + hull.equations[:, 2] <= 1e-8, axis=1)
        Z = np.full(len(tgt), np.nan); Z[mask] = LinearNDInterpolator(pts, allv)(tgt[mask])
        im = ax.pcolormesh(X, Y, Z.reshape(X.shape), shading='auto', cmap=cmap)
        closed = np.vstack([cell, cell[:1]]); ax.plot(closed[:, 0], closed[:, 1], 'k', lw=1)
        if show_kpoints: ax.scatter(pts[inside, 0], pts[inside, 1], s=3, c='w', alpha=.5)
        if path is not None:
            verts = np.asarray(path.kpoints if hasattr(path, 'kpoints') else [p[0] for p in path], float)
            labs = list(path.klabels if hasattr(path, 'klabels') else [p[1] for p in path])
            pp = verts @ bz.rlat*2*np.pi/BOHR2ANG @ frame
            ax.plot(pp[:, 0], pp[:, 1], 'w-', lw=1)
            for l, xy in zip(labs, pp):
                ax.annotate(l, xy, color='k', xytext=(4, 4), textcoords='offset points', annotation_clip=False)
        ax.set_xlim(xmin*1.06, xmax*1.06); ax.set_ylim(ymin*1.06, ymax*1.06)
        ax.set_aspect('equal'); ax.set_xlabel(r'$k_1$ (1/Å)'); ax.set_ylabel(r'$k_2$ (1/Å)')
        lab = (r'$\sum|\tilde d(\eta)|^2$' if quantity == 'weight' else r'$\sqrt{\sum|\tilde d(\eta)|^2}$')
        ax.set_title(rf'{lab},  $\eta$={eta}')
        fig.colorbar(im, ax=ax)
        if savefig: fig.savefig(savefig, bbox_inches='tight', dpi=200)
        return ax

    # ------------------------------------------------- Berry polarization / chi(n)
    def setup_berry(self, bz=None, overlap=None, field_direction=None):
        """
        Prepare the Berry-phase reconstruction of the polarization.
        bz      : FloquetBZ of the SAVE used in the RT run (default from save=...)
        overlap : path of ndb.Overlap (default: <folder>/<calc>/ndb.Overlap)
        field_direction : used to check that the BZ expansion respects the field
        """
        bz = self._bz(bz)
        if field_direction is not None: bz.check_field_symmetries(field_direction)
        if overlap is None:
            if self.samples.path is None: raise FloquetError('Give the path of ndb.Overlap')
            overlap = os.path.join(self.samples.path, 'ndb.Overlap')
        self.DIP_S = load_overlaps(overlap, bz.nk_full, self.samples.band_hi)
        if self.active_bands[-1] != bz.nocc:
            raise FloquetError('The propagated bands must include the last occupied band')
        return self.DIP_S

    def berry_polarization(self, C):
        """Berry-phase polarization [a.u., Cartesian] for coefficients C[t, k_IBZ, basis, state]."""
        if not hasattr(self, 'DIP_S'): self.setup_berry()
        bz = self.bz
        return berry_polarization(np.asarray(C)[:, bz.bz_to_ibz], self.DIP_S, bz.axes, bz.lattice_bohr,
                                  bz.volume, bz.spin_factor, bz.nocc, self.samples.band_lo, self.samples.band_hi)

    def polarization_from_floquet(self, tau):
        """Polarization at t0+tau computed from the Floquet-reconstructed states."""
        return self.berry_polarization(self.reconstruct_states(tau))

    def validate_berry(self, nldb, tol=1e-5, floor=1e-10):
        """
        Compare the Berry polarization of the stored RT states with the polarization
        written by yambo_nl (ndb.Nonlinear) at the same times. Returns metrics.
        """
        Praw = self.berry_polarization(self.C)
        Pref = stored_polarization_at(nldb, self.samples.times_au)
        dP = Praw - Praw[0]; dR = Pref - Pref[0]; rn = np.linalg.norm(dR, axis=0)
        if rn.max() <= 0: raise FloquetError('No induced polarization in ndb.Nonlinear')
        act = rn > floor*rn.max()
        err = np.linalg.norm(dP-dR, axis=0)/np.maximum(rn, 1e-300)
        m = {'active_components': [c for c, y in zip('xyz', act) if y],
             'relative_errors': {c: float(e) for c, e, y in zip('xyz', err, act) if y},
             'combined_error': _relnorm(dP[:, act]-dR[:, act], dR[:, act])}
        m['passed'] = bool(m['combined_error'] <= tol and np.all(err[act] <= tol))
        if not m['passed']:
            msg = f'Berry validation failed: {m}'
            if self.strict: raise FloquetError(msg)
            print('[WARNING] '+msg)
        return m

    def harmonic_response(self, nldb, order=2, fit_order=None, npoints=None, symmetry='auto',
                          sampling_tol=5e-4, max_refinements=3, validate=True, berry_tol=1e-5):
        """
        Harmonic response chi^(order)(-n omega; omega,...,omega) obtained from the Berry
        polarization of the Floquet-reconstructed states, with the same normalisation as
        the real-time analysis (Xn_from_sine). Returns chi (3,) complex [x,y,z
        field-contracted components, units (cm/statV)^(order-1)] and an info dict.

        nldb     : YamboNLDB of the same calculation (field, stored polarization)
        symmetry : 'none', 'residual' (average over field-preserving symmetries), or 'auto'
                   ('residual' if yambo_nl symmetrized the polarization)
        The number of time points is doubled until chi changes less than sampling_tol.
        """
        from yambopy.nl.external_efield import Divide_by_the_Field
        from yambopy.units import SVCMm12VMm1, AU2VMm1
        if order > self.nf: raise FloquetError(f'order {order} exceeds FLOrder={self.nf}')
        if len(nldb.Efield) != 1: raise FloquetError('One frequency per calculation expected')
        field = nldb.Efield[0]
        if str(field['name']).strip().upper() not in ('SIN', 'SOFTSIN'):
            raise FloquetError('Only SIN/SOFTSIN fields are supported')
        if not np.isclose(float(field['freq_range'][0]), self.omega_au, rtol=1e-6):
            raise FloquetError('ndb.Nonlinear frequency differs from ndb.RT_V_bands')
        if not hasattr(self, 'DIP_S'): self.setup_berry(field_direction=field['versor'])
        info = {}
        if validate: info['berry'] = self.validate_berry(nldb, tol=berry_tol)
        if symmetry == 'auto':
            symmetry = 'residual' if self.samples.path and read_symmetrize_flag(self.samples.path) else 'none'
        proj = polar_projector(self.bz, field['versor'], symmetry)
        fit_order = max(4, 2*order) if fit_order is None else fit_order
        n = npoints or max(65, 4*self.nf*self.nmodes+1, 2*fit_order+1)
        n += (n % 2 == 0)
        unit = np.power(SVCMm12VMm1/AU2VMm1, order-1)
        norm = Divide_by_the_Field(field, order)*unit

        def sample(npt):
            tau = np.arange(npt)*self.period_au/npt
            P = self.polarization_from_floquet(tau)
            P = P[0] + (P-P[0]) @ proj.T
            coef = harmonic_fit(self.t0_au+tau, P-P[0], self.omega_au, fit_order, float(field['initial_time']))
            return tau, P, coef, coef[order]*norm

        tau, P, coef, chi = sample(n); changes = []
        for _ in range(max_refinements):
            n2 = 2*n+1; tau2, P2, coef2, chi2 = sample(n2)
            scale = max(float(abs(chi2).max()), 1e-300)
            ch = float(np.max(abs(chi2-chi))/scale); changes.append(ch)
            n, tau, P, coef, chi = n2, tau2, P2, coef2, chi2
            if ch <= sampling_tol: break
        else:
            print(f'[WARNING] chi({order}) not converged with the time sampling: changes {changes}')
        info.update(symmetry=symmetry, projector=proj, npoints=n, sampling_changes=changes,
                    fit_order=fit_order, times_au=self.t0_au+tau, polarization=P, raw_harmonics=coef,
                    energy_eV=self.omega_au*ha2ev)
        return chi, info


##################################################################################
### PART II : Berry phase, overlaps and harmonic fit (module functions)
##################################################################################

def load_overlaps(filename, nk_full, nb):
    """
    Overlaps <u_k,m | u_k+dk,n> from ndb.Overlap (variable DIP_S), returned as
    (nk_full, 6, nb, nb) complex; directions 0,1,2 (forward neighbours along b1,b2,b3)
    are the ones used for the Berry phase.
    """
    if not os.path.isfile(filename):
        raise FloquetError(f'Missing {filename}: the Berry reconstruction needs ndb.Overlap')
    with Dataset(filename) as db:
        if 'DIP_S' not in db.variables: raise FloquetError(f'{filename} has no DIP_S variable')
        var = db.variables['DIP_S']; a = np.asarray(var[:], float); sh = a.shape
    if len(sh) != 6: raise FloquetError(f'DIP_S must have 6 dimensions, found {sh}')
    if sh == (1, nk_full, 6, nb, nb, 2):
        order = [1, 2, 3, 4, 5, 0]
    else:   # fall back on sizes
        ic = [i for i, n in enumerate(sh) if n == 2]; ik = [i for i, n in enumerate(sh) if n == nk_full]
        isp = [i for i, n in enumerate(sh) if n == 1]
        if len(ic) != 1 or len(ik) != 1 or len(isp) != 1 or nb == 6:
            raise FloquetError(f'Cannot decode DIP_S of shape {sh}')
        rest = [i for i in range(6) if i not in (ic[0], ik[0], isp[0])]
        idir = [i for i in rest if sh[i] == 6]
        bands = [i for i in rest if i not in idir and sh[i] == nb]
        if len(idir) != 1 or len(bands) != 2: raise FloquetError(f'Cannot decode DIP_S of shape {sh}')
        order = [ik[0], idir[0], bands[0], bands[1], ic[0], isp[0]]
    x = a.transpose(order)[..., 0]
    return x[..., 0] + 1j*x[..., 1]


def _full_occupied(C, nocc, band_lo, band_hi):
    """Add the frozen occupied bands below NLBands (identity) to the propagated ones."""
    V = np.zeros(C.shape[:-2] + (band_hi, nocc), complex)
    for b in range(band_lo-1): V[..., b, b] = 1.0
    V[..., band_lo-1:band_hi, band_lo-1:nocc] = C
    return V


def berry_polarization(C, DIP, axes, lattice, volume, spin, nocc, band_lo, band_hi):
    """
    King-Smith-Vanderbilt polarization along the strings of the regular grid,
    unwrapped in time on every string. C is [t, k_fullBZ, basis, state];
    returns the Cartesian polarization (nt, 3) in atomic units.
    """
    C = np.asarray(C, complex); nt, nk = C.shape[:2]
    if DIP.shape[:2] != (nk, 6): raise FloquetError('Overlap and k-grid dimensions disagree')
    phases = {a: np.empty((nt, len(info['cycles']))) for a, info in axes.items() if info['N'] > 1}
    for it in range(nt):
        V = _full_occupied(C[it], nocc, band_lo, band_hi)
        for a, ph in phases.items():
            info = axes[a]; vp = V[info['neighbor']]
            link = np.einsum('kbi,kbc,kcj->kij', V.conj(), DIP[:, a].transpose(0, 2, 1), vp, optimize=True)
            sign, logabs = np.linalg.slogdet(link)
            if np.any(~np.isfinite(logabs)) or np.any(logabs < np.log(1e-14)):
                raise FloquetError(f'Near-singular Berry link at time {it}, direction {a}')
            cyc = np.asarray(info['cycles'], int)
            ph[it] = np.angle(np.prod(sign[cyc], axis=1))
    p_red = np.zeros((nt, 3))
    for a, ph in phases.items():
        p_red[:, a] = spin/(2*np.pi)*np.unwrap(ph, axis=0).mean(axis=1)
    return p_red @ np.asarray(lattice, float)/volume


def polar_projector(bz, field, setting='none'):
    """Average of the unitary symmetries that leave the field invariant (polar-vector projector)."""
    if setting == 'none': return np.eye(3)
    if setting != 'residual': raise ValueError("symmetry must be 'none', 'residual' or 'auto'")
    e = np.asarray(field, float); e = e/np.linalg.norm(e)
    keep = [s for s, anti in zip(bz.sym_car, bz.antiunitary) if not anti and np.allclose(s @ e, e, atol=2e-6)]
    if not keep: raise FloquetError('No unitary field-preserving symmetry found')
    P = np.mean(keep, axis=0)
    if not np.allclose(P @ P, P, atol=2e-6): raise FloquetError('Symmetry average is not a projector')
    return P


def harmonic_fit(times_au, signal, omega, order, initial_time=0.0):
    """Least-squares fit S(t) = X0 + sum_n [X_n e^{-i n w t} + c.c.]; returns X_0..X_order."""
    T = np.asarray(times_au, float) - float(initial_time)
    M = np.zeros((len(T), 2*order+1), complex); M[:, 0] = 1.0
    for n in range(1, order+1):
        M[:, n] = np.exp(-1j*n*omega*T); M[:, n+order] = np.exp(+1j*n*omega*T)
    return np.linalg.lstsq(M, np.asarray(signal, complex), rcond=1e-12)[0][:order+1]


def stored_polarization_at(nldb, times_au):
    """Polarization written by yambo_nl interpolated (cubic spline) at the given times."""
    from scipy.interpolate import CubicSpline
    t = np.asarray(nldb.IO_TIME_points, float).ravel()
    P = np.asarray(nldb.Polarization[0], float).T
    if P.shape != (len(t), 3): raise FloquetError(f'Unexpected polarization shape {P.shape}')
    norms = np.linalg.norm(P-P[0], axis=1)
    if norms.max() == 0 or (len(t) > 20 and np.count_nonzero(norms > norms.max()*1e-12) <= 2):
        raise FloquetError('Stored polarization is empty or sparse (restarted run?)')
    tt = np.asarray(times_au, float)
    if tt.min() < t.min()-1e-9 or tt.max() > t.max()+1e-9:
        raise FloquetError('V_bands records fall outside the stored polarization times')
    return CubicSpline(t, P, axis=0)(np.clip(tt, t.min(), t.max()))


def read_symmetrize_flag(path):
    """True if the yambo_nl report says that current and polarization were symmetrized."""
    reps = sorted(glob.glob(os.path.join(path, 'r-*nloptics*')), key=os.path.getmtime)
    if not reps: return False
    with open(reps[-1], errors='replace') as f: text = f.read()
    hit = re.findall(r'Symmetrize current and pol\s*:\s*(yes|no)', text)
    return bool(hit) and hit[-1] == 'yes'


def floquet_harmonic_spectrum(calcs, folder='.', save='SAVE', order=2, verbose=False, **kwargs):
    """
    chi^(order) at the frequencies of a set of single-frequency Floquet runs.
    calcs : list of calculation folders (e.g. sorted(glob('FLQ_*')))
    Returns energies (eV, sorted) and chi (nfreq, 3).
    """
    from yambopy.dbs.nldb import YamboNLDB
    bz = FloquetBZ(save); rows = []
    for c in calcs:
        fa = FloquetAnalysis(folder=folder, calc=c, save=save, verbose=verbose)
        nl = YamboNLDB(folder=folder, calc=c)
        fa.setup_berry(bz=bz, field_direction=nl.Efield[0]['versor'])
        chi, info = fa.harmonic_response(nl, order=order, **kwargs)
        rows.append((info['energy_eV'], chi))
        print(f'{c}: {info["energy_eV"]:.4f} eV  chi = {chi}')
    rows.sort(key=lambda r: r[0])
    return np.array([r[0] for r in rows]), np.array([r[1] for r in rows])


##################################################################################
### PART III : legacy single-state API (unchanged, kept for backward compatibility)
###            Note: the secant/ratio search of VbPP treats one propagated state
###            at a time; FloquetAnalysis (Part I) is the recommended interface.
##################################################################################
#
# This class postprocess data from the ndb.V_bands and the ndb.V_bands_K_section databases 
# 
# 
class VbPP():

    def __init__(self,vb_db,kpt,band,vb_pars=None): 
        """vb_db: VBdb object (from YamboVbandsDB class)
           freq : float, field frequency
           vb_pars: list of floats
        """
        self.n_timesteps = vb_db.n_timesteps -1
        self.tvecs = self.read_tvecs(vb_db,kpt,band)
        times_= []
        for i in range(self.n_timesteps):
            times_.append(vb_db.times[i]/fs2aut) 
        self.times = times_ 
        self.kpt   = kpt
        self.band  = band
        self.basis_index = vb_db.basis_index
        self.basis_size = vb_db.basis_size
        self.ks_ev = get_bands(kpt=kpt,band=band)
        freq =  get_frequency(vb_db.vb_path)
        self.freq  =freq*ha2ev 
        self.period = 2.0*np.pi/freq/fs2aut #in fs
        self.qe_thrs,self.err_thrs,self.step,self.max_iter,self.max_fl_mode = self.set_conv_parameters(vb_db,vb_pars=vb_pars)
        self.tot_fl_modes = self.max_fl_mode * 2 + 1
        self.exp_mat_long,self.exp_mat_m1 = self.get_exp_matrix()

    def __str__(self):
        """
        Print all info of the class
        """
        s="\n * * * Vbands PP class  * * * \n\n"
        s+="Selected Kpt  : "+str(self.kpt)+"\n"
        s+="Selected Band : "+str(self.band)+"\n"
        s+="KS energy     : "+str(self.ks_ev)+" [eV]\n"
        s+="Basis index   : "+str(self.basis_index)+"\n"
        s+="\n =================================\n\n"
        s+="Field freq    : "+str(self.freq)+" [eV] \n"
        s+="Field period  : "+str(self.period)+" [fs] \n"
        s+="\n =================================\n\n"
        s+="QEnergy thrsh : "+str(self.qe_thrs)+"\n"
        s+="Error thrsh   : "+str(self.err_thrs)+"\n"
        s+="Conv step     : "+str(self.step)+"\n"
        s+="Max iteration : "+str(self.max_iter)+"\n"
        s+="\n =================================\n\n"
        s+="Max Floq mode : "+str(self.max_fl_mode)+"\n"
        s+="Tot Floq modes: "+str(self.tot_fl_modes)+"\n"
        s+="\n =================================\n\n"
        for i,v in enumerate(self.tvecs):
            s+="Time: "+str(self.times[i])+" fs:\n"
            s+='bnd   c.real             c.imag\n'
            for j in range(len(v)):
                idx = self.basis_index[0] + j
                s+=str(idx)+"   "+str(v[j].real)+"   "+str(v[j].imag)+"\n"
        return s

    def read_tvecs(self,vb_db,kpt,band):
        list_of_evecs = []
        for it in range(vb_db.n_timesteps-1):
            evec=np.zeros(vb_db.basis_size,dtype=complex)
            
            for i in range(vb_db.basis_size):
                evec[i] = vb_db.tvecs[it][kpt-1,band-1,i]
            list_of_evecs.append(evec)
        return list_of_evecs

    def set_conv_parameters(self,vb_db,vb_pars=None): # would be nice to be a dictionary
        if vb_pars is None:
            qe_thrs=1e-7
            err_thrs=1e-8
            step=0.01
            max_iter=300
            max_fl_mode = vb_db.fl_order
        else:
            qe_thrs,err_thrs,step,max_iter,max_fl_mode = vb_pars
            max_fl_mode = min(max_fl_mode,vb_db.fl_order)
        list_pars = [qe_thrs,err_thrs,step,max_iter,max_fl_mode]
        
        return list_pars
        
    def get_exp_matrix(self):
        """uses variables of class to set up a call to the
           external build_exp_matrix
        """
        mat, mat_m1 = build_exp_matrix(
                 listof_times=self.times,
                 max_fl_mode=self.max_fl_mode,
                 freq=self.freq)
        return mat, mat_m1

# Processing part
    
    def calc_pvecs(self,qe_ev):
        """Returns the periodic part of the
           Floquet basis functions
        """
        mat_of_pvecs = np.zeros((len(self.tvecs),self.basis_size),dtype=complex)
        for i,v in enumerate(self.tvecs):
            mat_of_pvecs[i,:] = np.exp(+1j * qe_ev * self.times[i] / hbar_eVfs) * v

        return mat_of_pvecs


    def calc_fvecs(self,qe_ev,mat_of_pvecs=None):
        """this will calculate the Floquet vectors for a given qe_evx
           and return """
        if self.exp_mat_m1 is None:
            raise AttributeError("You need to initialize the FL space first")
        if mat_of_pvecs is None:
            mat_of_pvecs = self.calc_pvecs(qe_ev)

        mat_of_fvecs = np.matmul(self.exp_mat_m1,mat_of_pvecs[:self.tot_fl_modes,:])

        return mat_of_fvecs

    def recalc_pvecs_via_fl(self,qe_ev=None,mat_of_pvecs=None,mat_of_fvecs=None):
        """Function to calculate pVecs via the obtained fVecs at all times steps,
           including those not used to generate those fVecs, i.e., outside the
           first period considered. This allow to determine whether the qe used
           truly makes the tVecs periodic
        """
        if qe_ev is None and mat_of_pvecs is None:
            raise ValueError("Provide either qe or pVecs")
        if qe_ev is None and mat_of_fvecs is None:
            raise ValueError("Provide either qe or fVecs")
        if mat_of_pvecs is None:
            mat_of_pvecs = self.calc_pvecs(qe_ev)
        if mat_of_fvecs is None:
            mat_of_fvecs = self.calc_fvecs(qe_ev,mat_of_pvecs=mat_of_pvecs)

        mat_of_recalc_pvecs = np.matmul(self.exp_mat_long,mat_of_fvecs)

        return mat_of_recalc_pvecs


    def run_nl2fl(self,qe_ev=None,tag=None,iter_num=None): #tag to be reviewed
        """run #TODO
        """
        if qe_ev is None:
            qe_ev = self.ks_ev

        if tag is None:
            tag = datetime.today().strftime('%Y%m%d-%H.%M.%S')
            if iter_num is not None:
                tag += '_iter'+str(iter_num)

        fl_eigenvectors = FLeigenvectors(self.period,self.freq,self.max_fl_mode,self.times,tag)

        nl_in = self.calc_pvecs(qe_ev=qe_ev)
        fl_out = self.calc_fvecs(qe_ev,mat_of_pvecs=nl_in)
        nl_out = self.recalc_pvecs_via_fl(mat_of_pvecs=nl_in,mat_of_fvecs=fl_out)
        err = np.sum(np.abs(nl_out - nl_in))
        
        fl_eigenvectors.store_results(nl_in,nl_out,fl_out,qe_ev,err)

        return fl_eigenvectors

    def qe_from_ratio(self,qe_ref=None,component=None):
        """Floquet quasienergy directly from the ratio of the raw coefficients
           at times differing by one period.

           For a Floquet state c_j(t) = exp(-i*qe*t/hbar) * p_j(t) with p_j(t+T)=p_j(t),
           hence  c_j(t+T)/c_j(t) = exp(-i*qe*T/hbar)   for every j and t.
           All (t,j) pairs are combined in the |c|^2-weighted estimator
               z = sum_{t,j} conj(c_j(t)) * c_j(t+T)   ->   qe = -hbar*arg(z)/T
           (so tiny, noisy components do not spoil the result).
           The phase is defined only modulo 2*pi, i.e. qe is defined modulo hbar*omega:
           the branch closest to qe_ref (default: KS energy) is returned.

           component : None (use all basis states) or 0-based index in the basis,
                       to use the ratio of a single coefficient only.
           Returns (qe_ev, info) with info a dict of diagnostics:
               eps_t     : qe estimated at each starting time (should be flat)
               spread    : std of eps_t [eV]
               coherence : |z| / sum|c(t)||c(t+T)|, equal to 1 for a perfect Floquet state
        """
        if qe_ref is None:
            qe_ref = self.ks_ev
        t = np.asarray(self.times)
        c = np.asarray(self.tvecs)
        T = self.period
        dt = np.diff(t)
        if t[-1] - t[0] < T:
            raise ValueError("Simulation shorter than one field period: cannot use the ratio method")

        shift = T/dt[0]
        if np.allclose(dt,dt[0],rtol=1e-6) and abs(shift-round(shift)) < 1e-3:
            # period is an integer number of time steps: no interpolation needed
            s = int(round(shift))
            c0, c1 = c[:len(t)-s], c[s:]
        else:
            print("\n\n WARNING!!! Period is not an integer number of the time steps \n\n")
            print("Time step :",str(dt[0])," fs ")
            print("Period    :",str(T)," fs ")
            # period not commensurate with the time grid: interpolate c(t+T)
            from scipy.interpolate import CubicSpline
            mask = (t + T) <= (t[-1] + 1e-9)
            c0 = c[mask]
            c1 = CubicSpline(t,c,axis=0)(t[mask]+T)

        if component is not None:
            c0, c1 = c0[:,component:component+1], c1[:,component:component+1]

        z_t = np.sum(np.conj(c0)*c1,axis=1)
        z   = np.sum(z_t)
        # phase is defined modulo 2pi -> qe modulo hbar*omega: choose branch closest to reference
        def _branch(eps0):
            return eps0 + np.round((qe_ref-eps0)/self.freq)*self.freq
        qe   = _branch(-hbar_eVfs*np.angle(z)/T)
        eps_t = _branch(-hbar_eVfs*np.angle(z_t)/T)
        coherence = np.abs(z)/np.sum(np.abs(c0)*np.abs(c1))
        info = {'eps_t':eps_t,'spread':np.std(eps_t),'coherence':coherence,
                'n_pairs':len(z_t)}
        return qe, info

    def plot_raw_tvecs(self,band_to_plot=1,outdir=None):
        """Plot the raw coefficients c_j(t) of the propagated state on KS state j
           (before removing the Floquet phase factor exp(-i qe t/hbar)).
           band_to_plot is 1-based position in the basis, as in plot_realtime.
        """
        c = np.array(self.tvecs)[:,band_to_plot-1]
        t = np.array(self.times)
        outdir = outdir or f'figs-raw_k{self.kpt}_b{self.band}'
        os.makedirs(outdir,exist_ok=True)
        fig,axes = plt.subplots(3,sharex=True)
        fig.set_size_inches(8.3,7.5)
        fig.suptitle(f'Raw coefficient on KS state {band_to_plot} (kpt {self.kpt}, band {self.band})')
        axes[0].plot(t,c.real,marker='.',color='tab:blue'); axes[0].set_ylabel(f'Re[ c_{band_to_plot} ]')
        axes[1].plot(t,c.imag,marker='.',color='tab:blue'); axes[1].set_ylabel(f'Im[ c_{band_to_plot} ]')
        axes[2].plot(t,np.abs(c),marker='.',color='tab:blue'); axes[2].set_ylabel(f'|c_{band_to_plot}|')
        axes[2].set_xlabel('Time (fs)')
        plt.show()
        #plt.savefig(f'{outdir}/fig-raw_coefficient_KS_state_{band_to_plot}.pdf')
        plt.close()

    def find_qe(self,qe_ev=None,tag=None,method='optimize',component=None):
        """method = 'optimize' : secant solver over run_nl2fl, minimising the
                                 error between nl_in and nl_out (original behaviour)
                    'ratio'    : quasienergy from c(t+T)/c(t) (see qe_from_ratio),
                                 no iteration; run_nl2fl is called once to get the
                                 Floquet coefficients and the periodicity error
           qe_ev     : initial guess ('optimize') / branch reference ('ratio')
        """
        if qe_ev   is None:
            qe_ev = self.ks_ev
        if method == 'ratio':
            return self._find_qe_ratio(qe_ev,tag=tag,component=component)
        if method != 'optimize':
            raise ValueError("method must be 'optimize' or 'ratio'")
        return self._find_qe_optimize(qe_ev,tag=tag)

    def _find_qe_ratio(self,qe_ev,tag=None,component=None):
        qe, info = self.qe_from_ratio(qe_ref=qe_ev,component=component)
        evecs = self.run_nl2fl(qe,tag=tag,iter_num=0)
        evecs.nr_it  = 0
        evecs.nr_acc = info['spread']
        evecs.ratio_info = info
        return evecs

    def _find_qe_optimize(self,qe_ev,tag=None):
        """ Secant solver to iterate over executions of run_nl2fl
            and minimize the error between nl_in and nl_out
        """
        lof_qe  = []
        lof_err = []

        # Iteration with guess
        evecs = self.run_nl2fl(qe_ev,tag=tag,iter_num=0)
        lof_qe.append(evecs.FL_qe)
        lof_err.append(evecs.err)
        if evecs.err < self.qe_thrs:
            evecs.nr_it = 0
            evecs.nr_acc = self.qe_thrs
            return evecs
        # Iteration with step
        evecs = self.run_nl2fl(qe_ev-self.step,tag=tag,iter_num=1)
        lof_qe.append(evecs.FL_qe)
        lof_err.append(evecs.err)
        _delta_qe = abs(lof_qe[-1]-lof_qe[-2])

        # Loop
        iter_num = 1
        while ((_delta_qe > self.qe_thrs) or (lof_err[-1] > self.err_thrs)) and (iter_num <= self.max_iter):
            iter_num += 1
            _derivative = (lof_err[-1]-lof_err[-2])/(lof_qe[-1]-lof_qe[-2])
            _qe = lof_qe[-1] - lof_err[-1]/_derivative
            evecs = self.run_nl2fl(_qe,tag=tag,iter_num=iter_num)
            lof_qe.append(evecs.FL_qe)
            lof_err.append(evecs.err)
            _delta_qe = abs(lof_qe[-1]-lof_qe[-2])

        evecs.nr_it = iter_num
        evecs.nr_acc = abs(lof_qe[-1]-lof_qe[-2])

        return evecs
#
# This class handles Floquet eigenvectors
# 
# 
class FLeigenvectors:

    def __init__(self,T,f,max_eta,listof_times,tag):
        self.period = T
        self.freq = f
        self.max_fl_mode = max_eta
        self.tot_fl_modes = 2 * max_eta + 1
        self.times = np.array(listof_times)
        self.tag = tag
        self.dir = 'figs-'+str(self.tag)

    @classmethod
    def FromArrayOnly(cls,array,tag):
        cls.FL_vecs = array
        cls.dir = 'figs-'+str(tag)
        cls.tot_fl_modes = array.shape[0]
        cls.max_fl_mode = int((cls.tot_fl_modes-1)/2)
        return cls

    def store_results(self,NL_in,NL_out,FL_vecs,FL_qe,err):
        self.NL_in = NL_in
        self.NL_out = NL_out
        self.FL_vecs = FL_vecs
        self.FL_qe = FL_qe
        self.err = err

    def calc_all_times_pVecs(self,t_step=0.0025,n_steps=669):
        M,alltimes = build_exp_matrix(
                     t0_fs=self.times[0],
                     t_step=t_step,
                     n_steps=n_steps,
                     max_fl_mode=self.max_fl_mode,
                     freq=self.freq,
                     l_inv=False)

        NL_out_alltimes = np.matmul(M,self.FL_vecs)
        return alltimes,NL_out_alltimes

    def plot_realtime(self,band_to_plot=1,t_step=0.0025):
        n_steps = int(self.period*4/t_step)+10
        X,Y=self.calc_all_times_pVecs(t_step=t_step,n_steps=n_steps)

        fig,axes=plt.subplots(2)
        fig.set_size_inches(8.3,5.8)
        fig.suptitle(f'Time-dependent projection over Kohn-Sham state: {band_to_plot+1}')

        axes[0].plot(X,Y[:,band_to_plot].real,label='FL_calculated',color='tab:blue')
        axes[0].plot(self.times,self.NL_in[:,band_to_plot].real ,label='NL_in' ,marker='o',linestyle='',color='tab:olive',ms=12)
        axes[0].plot(self.times,self.NL_out[:,band_to_plot].real,label='NL_out',marker='.',linestyle='',color='tab:blue',ms=11)
        axes[1].plot(X,Y[:,band_to_plot].imag,label='FL_calculated',color='tab:blue')
        axes[1].plot(self.times,self.NL_in[:,band_to_plot].imag ,label='NL_in' ,marker='o',linestyle='',color='tab:olive',ms=12)
        axes[1].plot(self.times,self.NL_out[:,band_to_plot].imag,label='NL_out',marker='.',linestyle='',color='tab:blue',ms=11)
        axes[1].set_xlabel('Time (fs)')
        axes[0].set_ylabel(f'Re[ d_{band_to_plot+1} ]')
        axes[1].set_ylabel(f'Im[ d_{band_to_plot+1} ]')
        os.system(f'if [ ! -d {self.dir} ]; then mkdir {self.dir};fi')
        plt.legend()
        plt.show()
        #plt.savefig(f'{self.dir}/fig-real_time_projection_over_KS_state_{band_to_plot+1}.pdf')
        #plt.close()

    def plot_realtime_raw(self,band_to_plot=1,t_step=0.0025):
        """Same as plot_realtime but for the coefficient BEFORE removing the Floquet
           phase, c(t) = exp(-i qe t/hbar) * p(t). Data points: raw NL coefficients;
           line: Floquet reconstruction exp(-i qe t/hbar) * sum_eta f_eta exp(-i eta w t).
        """
        n_steps = int(self.period*4/t_step)+10
        X,Y = self.calc_all_times_pVecs(t_step=t_step,n_steps=n_steps)
        phase_X = np.exp(-1j*self.FL_qe*X/hbar_eVfs)
        phase_t = np.exp(-1j*self.FL_qe*self.times/hbar_eVfs)
        C_fl  = Y[:,band_to_plot]*phase_X
        C_in  = self.NL_in[:,band_to_plot]*phase_t     # = raw tvecs
        C_out = self.NL_out[:,band_to_plot]*phase_t

        fig,axes = plt.subplots(3,sharex=True)
        fig.set_size_inches(8.3,7.5)
        fig.suptitle(f'Raw coefficient (with Floquet phase) on Kohn-Sham state: {band_to_plot+1}')
        for ax,f in zip(axes,[np.real,np.imag,np.abs]):
            ax.plot(X,f(C_fl),label='FL_calculated',color='tab:blue')
            ax.plot(self.times,f(C_in) ,label='NL_in (raw data)',marker='o',linestyle='',color='tab:olive',ms=12)
            ax.plot(self.times,f(C_out),label='NL_out',marker='.',linestyle='',color='tab:blue',ms=11)
        axes[0].set_ylabel(f'Re[ c_{band_to_plot+1} ]')
        axes[1].set_ylabel(f'Im[ c_{band_to_plot+1} ]')
        axes[2].set_ylabel(f'|c_{band_to_plot+1}|')
        axes[2].set_xlabel('Time (fs)')
        axes[0].legend()
        plt.show()
 #       os.makedirs(self.dir,exist_ok=True)
#        plt.savefig(f'{self.dir}/fig-real_time_raw_coefficient_KS_state_{band_to_plot+1}.pdf')
        plt.close()

    def plot_floquet(self,labels='+1'):
        _fl_vec = self.FL_vecs.reshape((1,self.FL_vecs.size),order='F')
        fig = plt.figure()
        fig.set_size_inches(8.3,3.9)
        ax = fig.gca()
        fig.suptitle('Projection over Floquet-Kohn-Sham space')

        cax = ax.matshow(abs(_fl_vec), cmap='YlOrBr', norm=pltcol.LogNorm(vmin=1.e-10,vmax=1.))
        fig.colorbar(cax,orientation='horizontal')

        lof_labels=[]
        for band in range(self.FL_vecs.shape[1]):
            for i in range(self.FL_vecs.shape[0]):
                tic='        '+str(i-self.max_fl_mode)
                lof_labels.append(tic)
        ax.set_xticks([x-0.5 for x in range(_fl_vec.shape[1])])
        ax.set_xticklabels(lof_labels)
        ax.set_yticks([])
        ax.set_yticklabels([])

        ax.text(0,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode+1,0]):.1e}', fontsize=15)
        ax.text(4,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode,0]):.1e}', fontsize=15)
        ax.text(10,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode-1,1]):.1e}', fontsize=15)
        ax.text(13.5,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode+1,1]):.1e}', fontsize=15)
        if (labels == '+2'):
            ax.text(16.5,-2.3, f'{abs(self.FL_vecs[self.max_fl_mode+2,1]):.1e}', fontsize=15)
        plt.grid(which='major',lw=1.,color='black')

        os.system(f'if [ ! -d {self.dir} ]; then mkdir {self.dir};fi')
        plt.savefig(f'{self.dir}/fig-FKS_projection.pdf')
        plt.close()

    def output(self):
        os.system(f'if [ ! -d {self.dir} ]; then mkdir {self.dir};fi')
        with open(f'{self.dir}/output-{self.tag}.dat','w') as f:
            for k,v in self.__dict__.items():
                f.write('--------------------------\n')
                f.write(f'Attribute: {k}\n')
                f.write('Values:\n')
                f.write(str(v)+'\n')

    def output_for_fortran(self,band,kpt,NL_band_1,file='fortran_input.txt'):
        if kpt == 1:
            filemode='w'
        else:
            filemode='a'
        with open(f'{file}',filemode) as f:
            for state in range(self.FL_vecs.shape[1]):
                for mode in range(self.FL_vecs.shape[0]):
                    _state = state + NL_band_1
                    f.write(f'FL_V_bands({_state},{mode+1},{band},{kpt},1) = {self.FL_vecs[mode,state].real}_SP + cI * {self.FL_vecs[mode,state].imag}_SP\n')
            f.write(f'FL_QE({band},{kpt},1) = {self.FL_qe/ha2ev}_SP\n')

##################################################################################
### AUXILIARY FUNCTIONS ###
##################################################################################

def build_exp_matrix(listof_times=None,t0_fs=None,t_step=None,n_steps=None,max_fl_mode=None,freq=None,l_inv=True):
  """builds matrix of exponentials with times
     and size given by self.times and self.tot_fl_modes
  """
  if max_fl_mode is None:
    raise ValueError("Missing total number of FL modes on input to build_exp_matrix")
  if listof_times is None and (t0_fs is None or t_step is None or n_steps is None):
    raise ValueError("Missing times on input to build_exp_matrix")

  if listof_times is None:
    listof_times = []
    for step in range(-int(n_steps*0.06),n_steps):
      t = t0_fs + t_step * step
      listof_times.append(t)

  tot_fl_modes = 2 * max_fl_mode + 1
  matrix = np.zeros((len(listof_times),tot_fl_modes),dtype=complex)
  for i,t in enumerate(listof_times):  # rows - time
   for j in range(tot_fl_modes): # columns - eta
     eta = j - max_fl_mode
     matrix[i,j]=np.exp(-1j * eta * (freq/hbar_eVfs) * t)

  if l_inv:
    inverse = np.linalg.inv(matrix[:tot_fl_modes,:])
    return matrix, inverse
  else:
    arrof_time = np.array(listof_times)
    return matrix, arrof_time

def get_frequency(db_path):
   ds=Dataset(db_path+'/ndb.Nonlinear')
   freq = float(ds['Field_Freq_1'][0])
   return freq

def get_bands(kpt=None,band=None):
  ds=Dataset('SAVE/ns.db1')
  KSevalues=ds['EIGENVALUES'][:]
  vb,vb_kpt,vbM=0,0,-111
  for _band in range(KSevalues.shape[2]):
    for _kpt in range(KSevalues.shape[1]):
      for _spin in range(KSevalues.shape[0]):
        if KSevalues[_spin,_kpt,_band] < 0. and KSevalues[_spin,_kpt,_band] > vbM:
           vb,vb__kpt,vbM = _band,_kpt,KSevalues[_spin,_kpt,_band]

  if kpt is None and band is None:
    return (KSevalues - vbM)*ha2ev
  else:
    return (KSevalues[0,kpt-1,band-1] - vbM)*ha2ev

##################################################################################
### FLOQUET BAND/PLOT FUNCTIONS ###
##################################################################################

def findallk_qe(vbdb,report_file=None,method='optimize'):
    from contextlib import nullcontext
    #
    ks_evk = np.zeros([vbdb.n_kpts,vbdb.n_vbands])
    fl_qek = np.zeros([vbdb.n_kpts,vbdb.n_vbands])
    fl_eig = np.zeros([vbdb.n_kpts,vbdb.n_vbands,vbdb.fl_order*2+1,vbdb.basis_size],dtype=complex) # COMPLEX CHANGE IT

    with open(report_file,'w') if report_file else nullcontext(sys.stdout) as f:
        for _bnd in range(vbdb.n_vbands):
            for _kpt in range(vbdb.n_kpts):
                vbs = VbPP(vbdb,_kpt+1,_bnd+1)
                opt_evecs = vbs.find_qe(method=method)
                ks_evk[_kpt,_bnd] = vbs.ks_ev
                fl_qek[_kpt,_bnd] = opt_evecs.FL_qe
                fl_eig[_kpt,_bnd,:,:] = opt_evecs.FL_vecs[:,:]
                f.write(f'Bnd {_bnd+1:2} @Kpt {_kpt+1:2}: KS eval [eV] = {ks_evk[_kpt,_bnd]:.6f} --- FL qe [eV] = {fl_qek[_kpt,_bnd]:.6f} -- Diff [eV] = {vbs.ks_ev-opt_evecs.FL_qe:.2e} -- It.: {opt_evecs.nr_it:2} -- err in period = {opt_evecs.err:.2e}\n')
    return fl_eig, fl_qek, ks_evk  

def calc_rho(vbdb,fl_eig,order,harmonic,kpt,bnd_1,bnd_2):
    '''
    vbdb     :: databases with info on vbands as output by YanboVbandsDB
    fl_eig   :: Floquet eigenvectors as output by findallk_qe
    order    :: integer, perturbation order of the field, Gamma (positive) 
    harmonic :: integer, harmonic order, gamma (non negative)
    kpt,band1,band2 :: integers, k point, and band indexes of the element of the density matrix

    RETURNS
    rho :: density matrix element (kpt,band1,band2) contributing to harmonic order gamma 
           for perturbation order Gamma, for all possible combinations of harmonic orders eta, nu of the fl_eig 
    '''
    # dimension ckecks
    args =locals()
    val =list(args.values())
    key = list(args.keys())
    for i in range(len(val)):
        if not(any([str(key[i]) == 'vbdb', str(key[i]) == 'fl_eig'])): 
            if not(isinstance(val[i], int)):
                raise ValueError("value for "+str(key[i])+" must be integer")        
            if val[i] <0:
                raise ValueError("value for "+str(key[i])+" must be nonnegative")
    if not((order - harmonic)%2 == 0 and ((order - harmonic)>=0)):
        raise ValueError("Harmonic order gamma not compatible with perturbation order Gamma: Gamma = gamma + 2n, where n is 0,1,...")
    if kpt >= vbdb.n_kpts:
        raise ValueError("Value for kpt must be smaller than "+str(vbdb.n_kpts))
    if any([bnd_1 < vbdb.basis_index[0], bnd_1 > vbdb.basis_index[1], bnd_2 < vbdb.basis_index[0], bnd_2 > vbdb.basis_index[1] ]):
        raise ValueError("Basis indexes must be in range "+ str( vbdb.basis_index))
    # basis index 
    b_1 = bnd_1 -  vbdb.basis_index[0]
    b_2 = bnd_2 -  vbdb.basis_index[0]
    # determine possible eta,nu values 
    if harmonic == order:
        eta = list(range(-harmonic,1))
    elif harmonic < order:
        eta = [-int((harmonic+order)/2),int((order-harmonic)/2)]
    nu = list(int(x + harmonic) for x in eta)
    c_dim = len(eta)
    rho = np.zeros(c_dim,dtype=complex)
    path = []
    for ii in range(c_dim):
         for n in range(vbdb.n_vbands):
             rho[ii] = rho[ii] + np.conj(fl_eig[kpt,n,eta[ii],b_2])*fl_eig[kpt,n,nu[ii],b_1]
         path.append([eta[ii],nu[ii]])
    return rho, path

def get_qebands_path(lat,path,fl_qek,ks_evk): # see how to put kwargs into it...
    from yambopy.plot.bandstructure import YambopyBandStructure
    from yambopy.lattice import red_car
    from yambopy.kpoints import get_path

    bands_kpoints, bands_indexes, path_car = get_path(lat.car_kpoints,lat.rlat,lat.sym_car,path,debug=False) 

    ks_bs = YambopyBandStructure(ks_evk[bands_indexes],bands_kpoints,kpath=path_car) 
    fl_bs = YambopyBandStructure(fl_qek[bands_indexes],bands_kpoints,kpath=path_car)

    return ks_bs, fl_bs
    
def get_qebands_interpolate(lat,path,fl_qek,ks_evk,lpratio=5,fermie=0,nelect = 0,verbose=1): 
    
    cell = (lat.lat, lat.red_atomic_positions, lat.atomic_numbers)
    trev_for_interp = lat.time_rev
    symrel = [sym for sym,trev in zip(lat.sym_rec_red,lat.time_rev_list) if trev==False ]
    kpoints = lat.red_kpoints
    _, _, path_car = get_path(lat.car_kpoints,lat.rlat,lat.sym_car,path,debug=False)
    band_kpoints_rlu = path.get_klist()[:,:3]
    band_kpoints = red_car(band_kpoints_rlu,lat.rlat) # why different from above???

    fl_eigens = np.zeros([2,np.shape(fl_qek)[0], np.shape(fl_qek)[1]])
    ks_eigens = np.zeros([2,np.shape(ks_evk)[0], np.shape(ks_evk)[1]])
    for ik in range(np.shape(fl_qek)[0]):
        fl_eigens[0,ik,:] = sorted(fl_qek[ik,:])
        ks_eigens[0,ik,:] = sorted(ks_evk[ik,:])
    skw_fl = SkwInterpolator(lpratio,kpoints,fl_eigens,fermie,nelect,cell,symrel,trev_for_interp,verbose=verbose)
    skw_ks = SkwInterpolator(lpratio,kpoints,ks_eigens,fermie,nelect,cell,symrel,trev_for_interp,verbose=verbose)

    fl_eigens_kpath = skw_fl.interp_kpts(band_kpoints_rlu).eigens[0]
    fl_bs = YambopyBandStructure(fl_eigens_kpath,band_kpoints,kpath=path_car)
    ks_eigens_kpath = skw_ks.interp_kpts(band_kpoints_rlu).eigens[0]
    ks_bs = YambopyBandStructure(ks_eigens_kpath,band_kpoints,kpath=path_car)

    return ks_bs, fl_bs

@add_fig_kwargs
def plot_2D_kdist(data,lat,nspin=-1,plt_cbar=False,shift_BZ=True,**kwargs):
    """
    2D scatterplot in the k-BZ of any real quantity which is a function of only the k-grid
    
    data:: real quantity function of k-grid
        - if plt_cbar colorbar is shown
        - if shift_BZ adjacent BZs are also plotted (default)
        - kwargs example: marker='H', s=300, cmap='viridis', etc.
        
    NB: THIS IS A COPY OF PLOT_DIPOLE except it is out of the dipole class and some minor detail
    """
    from yambopy.plot.plotting import add_fig_kwargs,BZ_Wigner_Seitz

    kpts =lat.car_kpoints
    rlat = lat.rlat

    # Input check
    if len(data)!=len(kpts):
        raise ValueError('Something wrong in data dimensions (%d data vs %d kpts)'%(len(data),len(kpts)))

    # Global plot stuff
    fig, ax = plt.subplots(1, 1)
    ax.add_patch(BZ_Wigner_Seitz(lat))

    if plt_cbar:
        if 'cmap' in kwargs.keys(): color_map = plt.get_cmap(kwargs['cmap'])
        else:                       color_map = plt.get_cmap('viridis')
    lim = 1.05*np.linalg.norm(rlat[0])
    ax.set_xlim(-lim,lim)
    ax.set_ylim(-lim,lim)

    # Reproduce plot also in adjacent BZs
    if shift_BZ:
        BZs = shifted_grids_2D(kpts,rlat)
        for kpts_s in BZs: plot=ax.scatter(kpts_s[:,0],kpts_s[:,1],c=data,**kwargs)
    else:
        plot=ax.scatter(kpts[:,0],kpts[:,1],c=data,**kwargs)

    if plt_cbar: cbar = fig.colorbar(plot)

    plt.gca().set_aspect('equal')

    return ax

    #if plt_show: plt.show()
    #else: print_string = "Plot ready.\nYou can customise adding savefig, title, labels, text, show, etc..."
