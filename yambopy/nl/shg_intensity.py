# Copyright (c) 2026, Myrta Gruening, Corey Nelson
# All rights reserved.
#
# This file is part of the yambopy project
# Calculate shg intensity from second-order susceptibility (yambo_nl + yambopy)
#

import os
import numpy as np
import yaml
from yambopy.units import AU2M, FREE_SPACE_PERM, speed_of_light_SI, hbar, electron_charge_SI, AU2KWCMm2

try:
    from refractiveindex import RefractiveIndexMaterial
except ImportError:
    RefractiveIndexMaterial = None



##########################################################################
# 1. Unit conversions
##########################################################################

HC_EV_M = 2.0 * np.pi * hbar * speed_of_light_SI / electron_charge_SI

# Gaussian(esu) -> SI conversion for a second-order susceptibility
ESU_TO_SI_CHI2 = 4.0 * np.pi / (speed_of_light_SI * 1e-4)


def omega_rad(omega_eV):
    "Photon energy in eV -> angular frequency in rad/s."
    return np.asarray(omega_eV, float) * electron_charge_SI / hbar


def chi2_supercell_to_sheet_SI(chi2_gaussian, Lz_SI):
    "Supercell-averaged chi^(2) (Gaussian) -> sheet chi^(2) in m^2/V."
    return np.asarray(chi2_gaussian) * ESU_TO_SI_CHI2 * Lz_SI


def sheet_to_bulk_chi2(chi2_sheet_SI, h_2D):
    "Effective bulk chi^(2) [m/V] from a sheet chi^(2) [m^2/V]."
    return np.asarray(chi2_sheet_SI) / h_2D


def intensity_au_to_SI(intensity_au):
    "Applied-field intensity in atomic units -> W/m^2."
    return float(intensity_au) * AU2KWCMm2 * 1e7


def nk_from_chi1_supercell(chi1_supercell, Lz_SI, h_2D):
    "Effective (n,k) of a 2D material from supercell chi^(1)."
    chi1_corrected = np.asarray(chi1_supercell) * Lz_SI / h_2D
    eps = 1.0 + 4.0 * np.pi * chi1_corrected
    eps1, eps2 = np.real(eps), np.imag(eps)
    mod = np.sqrt(eps1**2 + eps2**2)
    # np.maximum guards tiny negative arguments from floating-point noise
    # when eps2 ~ 0 (transparent regions), which would otherwise give NaN.
    n = np.sqrt(np.maximum(0.0, eps1 + mod) / 2.0)
    k = np.sqrt(np.maximum(0.0, -eps1 + mod) / 2.0)
    return n, k


def load_chi_order(path_or_folder, order, component="x"):
    """Load an o.YamboPy-X_probe_order_N file.
    Returns (omega_eV, chi) with chi complex (Gaussian units), per the
    yambopy output layout: E, Im(x), Re(x), Im(y), Re(y), Im(z), Re(z).
    component: "x" (default), "y", or "z".
    """
    if os.path.isdir(path_or_folder):
        path = os.path.join(path_or_folder, "o.YamboPy-X_probe_order_%d" % order)
    else:
        path = path_or_folder
    data = np.loadtxt(path)
    im = {"x": 1, "y": 3, "z": 5}[component]
    return data[:, 0], data[:, im + 1] + 1j * data[:, im]

def supercell_height_SI(lattice_db, axis=2):
    """Out-of-plane supercell height Lz in metres.
    axis selects the out-of-plane lattice vector (0=a, 1=b, 2=c);
    default 2=c. Assumes that lattice vector is aligned with a Cartesian
    axis, which is the usual setup; for a tilted cell the perpendicular
    spacing must be computed differently.
    """
    return lattice_db.lat[axis, axis] * AU2M

def field_intensity_SI(nonlinear_db_path, field_index=1):
    from netCDF4 import Dataset
    with Dataset(nonlinear_db_path, "r") as db:
        var = "Field_Intensity_%d" % field_index
        if var not in db.variables:
            raise KeyError("'%s' not found in %s; available: %s" % (var, nonlinear_db_path, list(db.variables)))
        intensity_au = float(db.variables[var][:])
    return intensity_au_to_SI(intensity_au)


class ChiLoader:
    
    COMPONENTS = ("x", "y", "z")
    
    def __init__(self, chidir, savedir=None, h_2D=0.0, lattice=None):
        self.chidir = chidir
        self.savedir = savedir
        self.h_2D = float(h_2D)
        
        self.omega_eV, self.order1 = self.load_order(1)
        _omega2, self.order2 = self.load_order(2)
        if not np.allclose(self.omega_eV, _omega2):
            raise ValueError("order 1 and order 2 files have different energy grids")
            
        if lattice is not None:
            self.lattice = lattice
        elif savedir is not None:
            from yambopy import YamboLatticeDB   # lazy import
            self.lattice = YamboLatticeDB.from_db_file(
                filename=os.path.join(savedir, "ns.db1"), Expand=False)
        else:
            self.lattice = None

        self.Lz = 0.0
        self.n = None
        self.k = None
        self.SHG_sheet = None
        self.SHG_bulk = None
    
    def load_order(self, order):
        path = os.path.join(self.chidir, "o.YamboPy-X_probe_order_%d" % order)
        data = np.loadtxt(path)
        # column pairs (Im, Re) per component: x=(1,2), y=(3,4), z=(5,6)
        chi = np.stack([data[:, i + 1] + 1j * data[:, i] for i in (1, 3, 5)], axis=1)
        return data[:, 0], chi

    def component(self, order, component="x"):
        "One Cartesian component of chi^(order) as a 1D complex array."
        idx = self.COMPONENTS.index(component)
        return (self.order1 if order == 1 else self.order2)[:, idx]

    def supercell_height_SI(self, axis=2):
        "Set and return the supercell height Lz in m."
        if self.lattice is None:
            raise ValueError("no lattice available: construct ChiLoader with savedir= or lattice=")
        self.Lz = supercell_height_SI(self.lattice, axis=axis)
        return self.Lz

    def nk_from_chi(self, component="x"):
        "Effective n, k of the monolayer from the run's own chi^(1)."
        if not (self.h_2D > 0 and self.Lz > 0):
            raise ValueError("nk_from_chi requires h_2D > 0 and Lz > 0 (set h_2D and call supercell_height_SI())")
        self.n, self.k = nk_from_chi1_supercell(self.component(1, component), self.Lz, self.h_2D)
        return self.n, self.k

    def chi2_supercell_to_sheet_SI(self):
        "Sheet chi^(2) in m^2/V for all components, shape (N, 3)."
        if not self.Lz > 0:
            raise ValueError("chi2_supercell_to_sheet_SI requires Lz > 0 (call supercell_height_SI() first)")
        self.SHG_sheet = chi2_supercell_to_sheet_SI(self.order2, self.Lz)
        return self.SHG_sheet

    def chi2_supercell_to_bulk_SI(self):
        "Bulk-equivalent chi^(2) in m/V (sheet / h_2D), shape (N, 3)."
        if not (self.h_2D > 0 and self.Lz > 0):
            raise ValueError("chi2_supercell_to_bulk_SI requires h_2D > 0 and Lz > 0")
        if self.SHG_sheet is None:
            self.chi2_supercell_to_sheet_SI()
        self.SHG_bulk = sheet_to_bulk_chi2(self.SHG_sheet, self.h_2D)
        return self.SHG_bulk

    def __str__(self):
        return ("\n * * * ChiLoader * * *\n"
                "Folder   : %s\nEnergies : %d (%.2f-%.2f eV)\n"
                "h_2D     : %g m\nLz       : %g m\n"
                % (self.chidir, len(self.omega_eV), self.omega_eV[0], self.omega_eV[-1], self.h_2D, self.Lz))
        


##########################################################################
# 2. refractiveindex.info database access
##########################################################################

DEFAULT_DB_ROOT = os.path.expanduser("~/.refractiveindex.info-database")

CATALOG = None
DB_ROOT = None


def load_catalog(db_root=None, force_reload=False):
    "Load (once) and return the refractiveindex.info catalog."
    global CATALOG, DB_ROOT
    root = db_root or os.environ.get("REFRACTIVEINDEX_DB", DEFAULT_DB_ROOT)
    if CATALOG is not None and not force_reload and root == DB_ROOT:
        return CATALOG
    catalog_path = os.path.join(root, "catalog-nk.yml")
    if not os.path.isfile(catalog_path):
        raise FileNotFoundError(
            "refractiveindex.info database not found at '%s'.\n"
            "It is auto-downloaded on the first RefractiveIndexMaterial(...) "
            "call (internet required), e.g.\n"
            "    from refractiveindex import RefractiveIndexMaterial\n"
            "    RefractiveIndexMaterial(shelf='main', book='Ag', page='Johnson')\n"
            "or set REFRACTIVEINDEX_DB / pass db_root=... to point at an "
            "existing copy." % root)
    with open(catalog_path) as f:
        CATALOG = yaml.safe_load(f)
    DB_ROOT = root
    return CATALOG


def database_version(db_root=None):
    "Database snapshot stamp from its .version file (None if absent)."
    root = db_root or DB_ROOT or os.environ.get("REFRACTIVEINDEX_DB",
                                                 DEFAULT_DB_ROOT)
    path = os.path.join(root, ".version")
    if os.path.isfile(path):
        with open(path) as f:
            return f.read().strip()
    return None


def search_database(keyword, db_root=None, exact_book=False):
    catalog = load_catalog(db_root)
    keyword = keyword.lower()
    results = []
    for shelf_entry in catalog:
        if 'SHELF' not in shelf_entry:      # skip top-level dividers
            continue
        shelf = shelf_entry['SHELF']
        for book_entry in shelf_entry.get('content', []):
            if 'BOOK' not in book_entry:    # skip dividers
                continue
            book = book_entry['BOOK']
            book_name = book_entry.get('name', book)
            if exact_book and book.lower() != keyword:
                continue
            for page_entry in book_entry.get('content', []):
                if 'PAGE' not in page_entry:  # skip dividers
                    continue
                page = page_entry['PAGE']
                source = page_entry.get('name', page)
                text = ("%s %s %s %s" % (book, book_name, page, source)).lower()
                if exact_book or keyword in text:
                    results.append({'shelf': shelf, 'book': book, 'page': page,
                                    'material': book_name, 'source': source})
    return results


def load_material(result):
    "Turn one search-result dict into a RefractiveIndexMaterial."
    if RefractiveIndexMaterial is None:
        raise ImportError("the 'refractiveindex' package is required: pip install refractiveindex")
    return RefractiveIndexMaterial(shelf=result['shelf'], book=result['book'],
                                   page=result['page'])


def print_search(keyword, db_root=None, exact_book=False):
    "Print matches with their energy range in eV and return the list."
    results = search_database(keyword, db_root=db_root, exact_book=exact_book)
    if not results:
        print("No matches for '%s'." % keyword)
        return results
    print("%d match(es) for '%s':\n" % (len(results), keyword))
    for i, r in enumerate(results):
        try:
            mat = load_material(r)
            lo, hi = mat.get_wl_range(unit='eV')
            lo, hi = min(lo, hi), max(lo, hi)
            erange = "%.2f - %.2f eV" % (lo, hi)
        except Exception:
            erange = "range unavailable"
        print("[%d] %s" % (i, r['material']))
        print("     shelf='%s', book='%s', page='%s'"
              % (r['shelf'], r['book'], r['page']))
        print("     %s" % r['source'])
        print("     energy range: %s\n" % erange)
    return results


def get_n(mat, value, unit='eV'):
    return mat.get_refractive_index(value, unit=unit)


def get_k(mat, value, unit='eV'):
    try:
        k = mat.get_extinction_coefficient(value, unit=unit)
    except Exception:
        return 0.0
    if k is None:
        return 0.0
    k = np.asarray(k, dtype=float)
    k = np.where(np.isnan(k), 0.0, k)
    return float(k) if k.ndim == 0 else k


def get_epsilon(mat, value, unit='eV'):
    n = get_n(mat, value, unit=unit)
    k = get_k(mat, value, unit=unit)
    return (np.asarray(n) + 1j * np.asarray(k))**2

##########################################################################
# 3. Material objects
##########################################################################

class Substrate(object):
    def __init__(self, name, record_index=0, source=None, db_root=None):
        self.name = name
        results = search_database(name, db_root=db_root)
        if not results:
            raise ValueError("No database match for '%s'" % name)
        if source is not None:
            hits = [r for r in results if source.lower() in r['source'].lower()]
            if not hits:
                raise ValueError(
                    "No '%s' record with source containing '%s'. Available: %s"
                    % (name, source, [r['source'] for r in results]))
            self.record = hits[0]
        else:
            if not (0 <= record_index < len(results)):
                raise ValueError(
                    "record_index %d out of range: '%s' has %d matches "
                    "(use print_search('%s') to list them)"
                    % (record_index, name, len(results), name))
            self.record = results[record_index]
        self.source = self.record['source']
        self.material = load_material(self.record)

    def wl_range_eV(self):
        lo, hi = self.material.get_wl_range(unit='eV')
        return (min(lo, hi), max(lo, hi))   # eV order flips vs wavelength

    def covers(self, energy_eV):
        e = np.asarray(energy_eV, float)
        lo, hi = self.wl_range_eV()
        return bool(np.all((e >= lo) & (e <= hi)))

    def check_range(self, energy_eV):
        if not self.covers(energy_eV):
            lo, hi = self.wl_range_eV()
            raise ValueError("%s (%s): requested energy outside dataset range %.2f-%.2f eV" % (self.name, self.source, lo, hi))

    def n(self, energy_eV):
        self.check_range(energy_eV)
        return get_n(self.material, energy_eV)

    def k(self, energy_eV):
        self.check_range(energy_eV)
        return get_k(self.material, energy_eV)

    def complex_index(self, energy_eV):
        return np.asarray(self.n(energy_eV)) + 1j * np.asarray(self.k(energy_eV))

    def epsilon(self, energy_eV):
        return self.complex_index(energy_eV)**2

    def __str__(self):
        lo, hi = self.wl_range_eV()
        return ("\n * * * Substrate * * *\n"
                "Material : %s\nSource   : %s\nValid    : %.2f-%.2f eV\n" % (self.name, self.source, lo, hi))


class SimulatedMaterial(object):

    def __init__(self, omega_eV, chi1_supercell, Lz_SI, h_2D,
                 name="2D material"):
        self.name = name
        self.h_2D = h_2D
        self.E = np.asarray(omega_eV, float)
        n, k = nk_from_chi1_supercell(np.asarray(chi1_supercell), Lz_SI, h_2D)
        self.n = np.asarray(n, float)
        self.k = np.asarray(k, float)
        self.lo, self.hi = self.E.min(), self.E.max()

    def n(self, e):
        return np.interp(np.asarray(e, float), self.E, self.n,
                         left=np.nan, right=np.nan)

    def k(self, e):
        return np.interp(np.asarray(e, float), self.E, self.k,
                         left=np.nan, right=np.nan)

    def complex_index(self, e):
        return np.asarray(self.n(e)) + 1j * np.asarray(self.k(e))

    def epsilon(self, e):
        return self.complex_index(e)**2

    def wl_range_eV(self):
        return (self.lo, self.hi)

    def covers(self, e):
        e = np.asarray(e, float)
        return bool(np.all((e >= self.lo) & (e <= self.hi)))

    def __str__(self):
        return ("\n * * * SimulatedMaterial * * *\n"
                "Material : %s\nSource   : chi(1)\nValid    : %.2f-%.2f eV\n"
                % (self.name, self._lo, self._hi))




##########################################################################
# 4. SHG intensity models
##########################################################################

class Stack(object):
    "Air / 2D / film(d) / substrate. structure-factor SHG model."


    def __init__(self, material_2D, film, substrate, film_thickness, h_2D):
        self.material_2D = material_2D
        self.film = film
        self.substrate = substrate
        self.d = film_thickness
        self.h_2D = h_2D
        self.n0 = 1.0            # air

    def usable_omega(self, omega_eV):
        "Energies where every layer has data at both w and 2w."
        w = np.asarray(omega_eV, float)
        usable = np.ones(len(w), bool)
        for mat in (self.material_2D, self.film, self.substrate):
            lo, hi = mat.wl_range_eV()
            usable &= (w >= lo) & (w <= hi) & (2*w >= lo) & (2*w <= hi)
        return usable

    def r_ij(self, ni, nj):
        "Fresnel reflection coefficient between two media."
        return (ni - nj)/(ni + nj)

    def Rs(self, wavelength, n1, n2):
        "Film+substrate reflection without the 2D layer, including film interference.  Song et al. Eq. (4)."
        w_tilde = 2*np.pi/wavelength
        r01 = self.r_ij(self.n0, n1)
        r12 = self.r_ij(n1, n2)
        phase = np.exp(2j*w_tilde*n1*self.d)
        return (r01 + r12*phase)/(1 + r01*r12*phase)

    def R_total(self, wavelength, n2D, n1, n2):
        "Total reflection of air/2D/film/substrate.  Song et al. Eq. (3)."
        w_tilde = 2*np.pi/wavelength
        eta = -1j*self.h_2D*w_tilde*(n2D**2 - 1)/2
        r = -eta/(1 + eta)
        t = 1/(1 + eta)
        Rs = self.Rs(wavelength, n1, n2)
        return r + (Rs*t**2)/(1 - Rs*r)

    def structure_factor(self, omega_eV):
        "beta = (1 + R_w)^2 (1 + R_2w), Song et al. Eq. (2)"
        w = np.asarray(omega_eV, float)
        usable = self.usable_omega(w)
        beta = np.full(len(w), np.nan, dtype=complex)
        if not usable.any():
            print("Stack: no omega usable by all materials at omega and 2*omega")
            return beta
        wv = w[usable]
        wavelength = HC_EV_M/wv
        R_w = self.R_total(wavelength, self.material_2D.complex_index(wv), self.film.complex_index(wv), self.substrate.complex_index(wv))
        R_2w = self.R_total(wavelength/2, self.material_2D.complex_index(2*wv), self.film.complex_index(2*wv), self.substrate.complex_index(2*wv))
        beta[usable] = (1 + R_w)**2 * (1 + R_2w)
        return beta

    def shg_intensity(self, omega_eV, chi2_sheet, I_incident):
        "SHG intensity, Song et al. Eq. (1) for chi2_sheet in m^2/V"
        w = np.asarray(omega_eV, float)
        wavelength = HC_EV_M/w
        beta = self.structure_factor(w)
        chi_term = np.abs(2*np.pi*np.asarray(chi2_sheet)/wavelength)**2
        return (1/(2*FREE_SPACE_PERM*speed_of_light_SI)) * np.abs(beta)**2 * chi_term * I_incident**2

    def __str__(self):
        return ("\n * * * Stack * * *\n"
                "2D material : %s\nFilm        : %s (d = %.1f nm)\n"
                "Substrate   : %s\nh_2D        : %.2f nm\n"
                % (getattr(self.material_2D, 'name', '?'),
                   getattr(self.film, 'name', '?'), self.d*1e9,
                   getattr(self.substrate, 'name', '?'), self.h_2D*1e9))


class WoodwardModel:
# used for sheet
  def __init__(self, material_2D, substrate, h_2D):
    self.material_2D = material_2D
    self.substrate = substrate
    self.h_2D = h_2D

  def omega_rad(self, omega_eV):
    return np.asarray(omega_eV, float) * electron_charge_SI / hbar

  def usable_omega(self, omega_eV):
     w = np.asarray(omega_eV, float)
     lo, hi = self.substrate.wl_range_eV()
     return (w >= lo) & (w <= hi) # Checks if the selected refractive index range covers the users data

  def sheet_intensity(self, omega_eV, chi2_sheet, I_incident):
    w = np.asarray(omega_eV, float)
    usable = self.usable_omega(w)
    I = np.full(len(w), np.nan, float)
    if not usable.any():
         print("WoodwardModel: no omega usable by the substrate."); return I
    wv = w[usable] # Valid omega values
    omega = self.omega_rad(wv)
    n_complex = np.asarray(self.substrate.complex_index(wv))
    if np.any(np.abs(np.imag(n_complex)) > 0.01):
      print("WARNING: Woodward sheet model assumes a TRANSPARENT substrate, but this one absorbs (k>0). Use the structure-factor Stack instead.") #
    n_q = np.real(n_complex)
    chi = np.asarray(chi2_sheet)[usable] if np.ndim(chi2_sheet) else chi2_sheet

    num = 32 * omega**2 * np.abs(chi)**2 * I_incident**2
    den = (n_q + 1)**6 * FREE_SPACE_PERM * speed_of_light_SI**3

    I[usable] = num / den
    return I

  def __str__(self):
        return (f"\n * * * WoodwardModel * * *\n"
                f"2D material : {getattr(self.material_2D,'name','?')}\n"
                f"Substrate   : {getattr(self.substrate,'name','?')}\n"
                f"h_2D        : {self.h_2D*1e9:.2f} nm\n")
    

class ClarkModel:
# used for the bulk
  def __init__(self, material_2D, substrate, h_2D):
    self.material_2D = material_2D
    self.substrate = substrate
    self.h_2D = h_2D


  def omega_rad(self, omega_eV):
    return np.asarray(omega_eV, float) * electron_charge_SI / hbar


  def usable_omega_bulk(self, omega_eV):
    w = np.asarray(omega_eV, float)
    lo, hi = self.material_2D.wl_range_eV()
    return (w >= lo) & (w <= hi) & (2*w >= lo) & (2*w <= hi) # Same as above with extra check for 2*omega, due to its use in the bulk model formula


  def bulk_intensity(self, omega_eV, chi2_bulk, I_incident): # Equation 4
    w = np.asarray(omega_eV, float)
    usable = self.usable_omega_bulk(w)
    I = np.full(len(w), np.nan, float)
    if not usable.any():
         print("ClarkModel: no omega usable by the 2D material at omega and 2*omega."); return I
    wv = w[usable]
    omega = self.omega_rad(wv)
    n_w = np.real(np.asarray(self.material_2D.complex_index(wv)))
    n_2w = np.real(np.asarray(self.material_2D.complex_index(2*wv)))
    chi = np.asarray(chi2_bulk)[usable] if np.ndim(chi2_bulk) else chi2_bulk
    num = omega**2 * np.abs(chi)**2 * self.h_2D**2 * I_incident**2
    den = 2 * n_w**2 * n_2w * FREE_SPACE_PERM * speed_of_light_SI**3
    I[usable] = num / den
    return I


  def __str__(self):
        return (f"\n * * * ClarkModel * * *\n"
                f"2D material : {getattr(self.material_2D,'name','?')}\n"
                f"Substrate   : {getattr(self.substrate,'name','?')}\n"
                f"h_2D        : {self.h_2D*1e9:.2f} nm\n")
    
    
    
