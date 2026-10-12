# Using `shg_intensity`: a practical tutorial

This tutorial shows how to use the `shg_intensity` module
(`yambopy/nl/shg_intensity.py`) to turn the output of a `yambo_nl` run into a
second-harmonic-generation (SHG) intensity spectrum you can compare with
experiment. It is organised around tasks: load a run, convert the numbers, set
up the materials, compute the intensity, check the result. Short explanations
of why each step exists are kept to a few lines; the code blocks are the
main content and every snippet is runnable as shown.

Conventions: the module works in $SI$ units (lengths in $m$, intensities in
$W m^{-2}$, sheet susceptibilities in $m^2 V^{-1}$), photon energies are
passed in $eV$, and the raw yambo files are in Gaussian units; the conversion
to SI is one of the steps below.

Every snippet assumes:

```python
import numpy as np
from yambopy import *
from yambopy.nl.shg_intensity import *
```

---

## 0. Quick reference

The whole pipeline in one block: copy, point it at your run, and go. Sections
1–5 explain each line; Section 6 walks through it with real output.

```python
# load the run (reads the chi1 and chi2 files and the lattice database)
CHI = ChiLoader('.', 'SAVE', h_2D=0.65e-9)

# convert
CHI.supercell_height_SI()          # -> CHI.Lz
CHI.chi2_supercell_to_sheet_SI()   # -> CHI.SHG_sheet
CHI.chi2_supercell_to_bulk_SI()    # -> CHI.SHG_bulk
CHI.nk_from_chi()                  # -> CHI.n, CHI.k

# materials
MoS2   = SimulatedMaterial(CHI.omega_eV, CHI.component(1, 'x'),
                           CHI.Lz, CHI.h_2D, name="MoS2")
silica = Substrate("SiO2", record_index=8)
si     = Substrate("Si",   record_index=200)
I0     = field_intensity_SI('SAVE/ndb.Nonlinear')

# intensity of MoS2 / SiO2(285 nm) / Si
stack = Stack(MoS2, silica, si, film_thickness=285e-9, h_2D=CHI.h_2D)
I_shg = stack.shg_intensity(CHI.omega_eV, CHI.SHG_sheet[:, 0], I0)

# cross-check (should be ~1.0)
stack0 = Stack(MoS2, silica, silica, film_thickness=0.0, h_2D=CHI.h_2D)
I_w    = WoodwardModel(MoS2, silica, CHI.h_2D).sheet_intensity(
             CHI.omega_eV, CHI.SHG_sheet[:, 0], I0)
print(np.nanmedian(stack0.shg_intensity(CHI.omega_eV, CHI.SHG_sheet[:, 0], I0) / I_w))
```

---

## 1. Load a run

You need three things from a `yambo_nl` run analysed with `Xn_from_sine`: the
folder containing `o.YamboPy-X_probe_order_1` and `_2`, the run's `SAVE`
folder, and a value for the effective monolayer thickness.

```python
CHI = ChiLoader('.', 'SAVE', h_2D=0.65e-9)
print(CHI)     # summary: number of energies, energy range, h_2D, Lz
```

One `ChiLoader` object represents one run. Constructing it reads both
response orders; the `order_1` file (the linear susceptibility $\chi^{(1)}$)
and the `order_2` file (the second-order susceptibility $\chi^{(2)}$, the SHG
quantity), and it opens the lattice database. Both are needed later:
$\chi^{(2)}$ is the source of the SHG signal, and $\chi^{(1)}$ is what the
monolayer's refractive index is built from. `h_2D` is the effective monolayer
thickness in m (0.65 nm is the conventional $MoS_2$ value).

After construction you have:

```python
CHI.omega_eV        # the photon-energy grid of the run, in eV (length N)
CHI.order1          # chi1 spectrum
CHI.order2          # chi2 spectrum
```

`order1` and `order2` are arrays of shape (N, 3): one row per photon
energy (N photon energies), one column per Cartesian direction (x, y, z),
matching the file layout
`E[eV], Im(x), Re(x), Im(y), Re(y), Im(z), Re(z)`. Each entry is a single
complex number, the loader combines each Im/Re column pair into the normal complex number format (Re + i Im). All three directions are loaded because $\chi^{(2)}$ is
in general a vector-valued (3D) object; you choose a direction only when you use
it:

```python
chi2_x = CHI.component(2, 'x')     # 1D complex array, length N
chi2_y = CHI.component(2, 'y')     # the second in-plane component
CHI.order2[:, 0]                   # same as CHI.component(2, 'x')
```

**Do this immediately after loading:** it costs one line and catches broken
runs:

```python
assert np.allclose(CHI.component(2, 'z'), 0.0)
```

A flat monolayer (point group $D_{3h}$) cannot have an out-of-plane
$\chi^{(2)}$; the z columns must be numerically zero. 

---

## 2. Convert to physical quantities

The raw file values are in Gaussian units and are averaged over the whole
simulation supercell, which is mostly vacuum, so their magnitude depends on
how much vacuum padding the run used. Four method calls produce the quantities
you actually work with. Each stores its result as an attribute on `CHI` (and
also returns it):

```python
CHI.supercell_height_SI()          # stores CHI.Lz, the supercell height in m
CHI.chi2_supercell_to_sheet_SI()   # stores CHI.SHG_sheet, shape (N, 3), m^2/V
CHI.chi2_supercell_to_bulk_SI()    # stores CHI.SHG_bulk,  shape (N, 3), m/V
CHI.nk_from_chi()                  # stores CHI.n and CHI.k, length-N real arrays
```

Call `supercell_height_SI` first; the other three need `Lz` and will raise a
clear error telling you so if it is not set. What each quantity is for:

* **`CHI.Lz`** is read from `SAVE/ns.db1`, i.e. the value the calculation
  actually used that is never typed by hand, so it cannot be mistyped.
* **`CHI.SHG_sheet`** is the sheet susceptibility: the raw $\chi^{(2)}$
  multiplied by `Lz` (undoing the vacuum dilution) and converted to SI. This is
  the vacuum-independent quantity and the one you feed to the intensity models.
* **`CHI.SHG_bulk`** is **CHI.SHG_sheet** / **$h_{2D}$**, the bulk-equivalent value
  that papers quote in pm/V. Use it only to compare with literature, and
  remember it depends on the thickness convention: a paper using $h = 0.61$ nm
  will quote a different pm/V for identical physics.
* **`CHI.n`, `CHI.k`** are the monolayer's effective optical constants,
  computed from the run's own $\chi^{(1)}$ so that the linear and nonlinear
  parts of the analysis are consistent.

---

## 3. Set up the materials

The intensity models need optical constants for every layer. These come from
two places, wrapped in two classes with the same interface.

**The monolayer** comes from the simulation itself:

```python
MoS2 = SimulatedMaterial(CHI.omega_eV, CHI.component(1, 'x'),
                         CHI.Lz, CHI.h_2D, name="MoS2")
print(MoS2)                      # summary with validity range
MoS2.n(1.9), MoS2.k(1.9)         # n and k at 1.9 eV
```

**The substrate layers** come from the refractiveindex.info database
(auto-downloaded on first use):

```python
print_search("SiO2")             # lists available records with energy ranges
silica = Substrate("SiO2", record_index=8)
si     = Substrate("Si",   record_index=200)
print(silica)                    # provenance: source paper and validity range
```

A material name alone is ambiguous as the database holds many measured records
per material, covering different energy ranges. `print_search` shows the
options so the choice is explicit; the printed `source` is provenance for your
methods section. Two practical notes: selecting by `source="..."` is stable
against database updates while `record_index` is not, and `database_version()`
records which snapshot you used. `Substrate` raises an error if you ask for a
value outside its measured range (extrapolating data is silently wrong);
`SimulatedMaterial` returns NaN outside its grid instead, so plots survive.

**The pump intensity** is read from the run's field database, so the predicted
intensity is consistent with the $\chi^{(2)}$ extracted from the same field:

```python
I0 = field_intensity_SI('SAVE/ndb.Nonlinear')     # W/m^2
```

---

## 4. Compute the SHG intensity


For the standard experimental geometry, monolayer on a dielectric film on a
substrate, build a `Stack` and call `shg_intensity`:

```python
stack = Stack(MoS2, silica, si, film_thickness=285e-9, h_2D=CHI.h_2D)
I_shg = stack.shg_intensity(CHI.omega_eV, CHI.SHG_sheet[:, 0], I0)   # W/m^2
```

Note the explicit `[:, 0]`: the models take one component of the (N, 3)
sheet susceptibility, and you choose which (here x). Making that choice at the
call site keeps the tensor-component decision visible in your analysis instead
of hidden in a loader default.

The dielectric film acts as an interferometer, causing incoming light to constructively and destructively interfere with itself. That physics is isolated in the structure factor $\beta$, which you can inspect on its own:

```python
beta = stack.structure_factor(CHI.omega_eV)
print(abs(beta)**2)          # |beta|^2 ranges from 0 to 64
```


**Expect NaN at the top of the spectrum.** Every energy needs optical data at
both $\omega$ and $2\omega$ for every layer, so your usable fundamental range
ends at half the maximum tabulated energy of your narrowest material. Energies
that fail the test come back NaN (with a printed warning), never silently
extrapolated. To see the usable window:

```python
ok = stack.usable_omega(CHI.omega_eV)     # boolean array, True where valid
print(CHI.omega_eV[ok].min(), "to", CHI.omega_eV[ok].max(), "eV usable")
```

---

## 5. Check the result before you believe it

The module ships three independently derived intensity models precisely so you
can check them against each other. The standard check: **run it once for
every new dataset**. If a 2D material sits directly on a transparent substrate
(zero film thickness), `Stack` and `WoodwardModel` describe the same
physics through completely independent code paths:

```python
stack0 = Stack(MoS2, silica, silica, film_thickness=0.0, h_2D=CHI.h_2D)
I_st0  = stack0.shg_intensity(CHI.omega_eV, CHI.SHG_sheet[:, 0], I0)
I_wood = WoodwardModel(MoS2, silica, CHI.h_2D).sheet_intensity(
             CHI.omega_eV, CHI.SHG_sheet[:, 0], I0)

print("median ratio:", np.nanmedian(I_st0 / I_wood))    # expect ~1.0
```

The third model is the bulk thin-slab reference (Clark's Eq. (4)), useful for
thickness-scaling comparisons. Note it takes the **bulk-equivalent**
$\chi^{(2)}$, not the sheet value:

```python
I_bulk = ClarkModel(MoS2, silica, CHI.h_2D).bulk_intensity(
             CHI.omega_eV, CHI.SHG_bulk[:, 0], I0)
```

`WoodwardModel` warns if you hand it an absorbing substrate as its formula
assumes transparency; use `Stack` for absorbing substrates like Si.

---

## 6. A worked example, start to finish

This section runs the whole pipeline on one concrete case: **monolayer $MoS_2$
on 285 nm $SiO_2$ on Si**, a fundamental grid spanning 0.50–2.97 eV, and shows
the output at each stage so you have something to check your own run against.
The outputs shown are from a real run; your own numbers will differ, and that
is fine. 

**Load and sanity-check.**

```python
CHI = ChiLoader('mos2_shg', 'mos2_shg/SAVE', h_2D=0.65e-9)
print(CHI)
assert np.allclose(CHI.component(2, 'z'), 0.0)     # D3h: no out-of-plane chi2
```

```text
ChiLoader: 100 energies over 0.50–2.97 eV, 3 components (x, y, z)
h_2D = 6.50e-10 m,  Lz = 2.1167e-09 m (Lz is the supercell height found using CHI.supercell_height_SI())
```

The z assertion passing is your first confirmation the run is physical, if it
raises, stop and check the simulation before going further.

**Convert, and read off a number you can compare with the literature.**

```python
CHI.supercell_height_SI()
CHI.chi2_supercell_to_sheet_SI()
CHI.chi2_supercell_to_bulk_SI()
CHI.nk_from_chi()
```



```python
# Peak sheet |chi2| across the spectrum (a quick magnitude sanity check)
i = np.nanargmax(np.abs(CHI.SHG_sheet[:, 0])) # finds index of peak
print("peak sheet |chi2|:", abs(CHI.SHG_sheet[i, 0]), "m^2/V at",
     CHI.omega_eV[i], "eV")
print("bulk-equiv |chi2| there:", abs(CHI.SHG_bulk[i, 0]) * 1e12, "pm/V")
```

The bulk-equivalent value in pm/V is the quantity to sanity-check: monolayer
$MoS_2 \chi^{(2)}$ is expected in the hundreds of pm/V, so a peak of this
order confirms the unit chain is right. Remember (Section 2) that any pm/V
comparison is only as fixed as the assumed monolayer thickness; the same
physics quoted with a different $h$ gives a different number.

**Build the materials and compute the spectrum.**

```python
MoS2   = SimulatedMaterial(CHI.omega_eV, CHI.component(1, 'x'),
                           CHI.Lz, CHI.h_2D, name="MoS2")
silica = Substrate("SiO2", record_index=8)
si     = Substrate("Si",   record_index=200)
I0     = field_intensity_SI('mos2_shg/SAVE/ndb.Nonlinear')

stack = Stack(MoS2, silica, si, film_thickness=285e-9, h_2D=CHI.h_2D)
I_shg = stack.shg_intensity(CHI.omega_eV, CHI.SHG_sheet[:, 0], I0)

ok = stack.usable_omega(CHI.omega_eV)
print(CHI.omega_eV[ok].min(), "to", CHI.omega_eV[ok].max(), "eV usable")
```

```text
0.50 to 1.47 eV usable   (40 of 100 points)
```

Note what the coverage rule (Section 4) does here: with a fundamental grid to
2.97 eV, only fundamentals whose $2\omega$ stays inside every layer's data
range survive; so 40 of the 100 points make it through, and the usable window
is the lower half of the grid (up to $\approx 1.47$ eV, i.e. $2\omega \approx
2.94$ eV). If this line prints an empty or tiny range for your run, that is the
coverage limit, not a bug; extend the $\chi$ grid or pick wider-range
substrate records.

**Plot the result.**

```python
import matplotlib.pyplot as plt
plt.plot(CHI.omega_eV[ok], I_shg[ok])
plt.xlabel("fundamental photon energy (eV)")
plt.ylabel(r"$I(2\omega)$ (W m$^{-2}$)")
plt.title("SHG intensity — MoS$_2$ / 285 nm SiO$_2$ / Si")
plt.show()
```

You should see a curve tracking the $MoS_2$ resonances, modulated by the
$SiO_2/Si$ interference through the structure factor.

**Confirm you believe it.** The last thing to do is the
independent cross-check:

```python
stack0 = Stack(MoS2, silica, silica, film_thickness=0.0, h_2D=CHI.h_2D)
ratio  = (stack0.shg_intensity(CHI.omega_eV, CHI.SHG_sheet[:, 0], I0)
          / WoodwardModel(MoS2, silica, CHI.h_2D).sheet_intensity(
                CHI.omega_eV, CHI.SHG_sheet[:, 0], I0))
print("d=0 Stack/Woodward median ratio:", np.nanmedian(ratio))
```

```text
d=0 Stack/Woodward median ratio: 0.9611
```

A ratio near 1.0 means the two independent code paths agree, here 0.96, within a few percent of unity, which is the pass
you want. That is the point at which the spectrum above is worth trusting.

---

## 7. Common tasks

**Plot the spectrum:**

```python
import matplotlib.pyplot as plt
plt.plot(CHI.omega_eV, I_shg)
plt.xlabel("photon energy (eV)"); plt.ylabel(r"$I(2\omega)$ (W m$^{-2}$)")
plt.show()
```

**Compare two runs (e.g. monolayer vs bulk).** Make one loader per run; each
object keeps its own `Lz`, spectra, and results, so nothing can cross-contaminate:

```python
mono = ChiLoader('mono_run', 'mono_run/SAVE', h_2D=0.65e-9)
bulk = ChiLoader('bulk_run', 'bulk_run/SAVE', h_2D=0.65e-9)
```

**Scan the film thickness** to find where the structure factor peaks:

```python
for d in [90e-9, 285e-9, 300e-9]:
    b = Stack(MoS2, silica, si, film_thickness=d, h_2D=CHI.h_2D).structure_factor(CHI.omega_eV)
    print(d, np.nanmax(abs(b)**2))
```

**Look at the y component** (the second in-plane direction):

```python
I_shg_y = stack.shg_intensity(CHI.omega_eV, CHI.SHG_sheet[:, 1], I0)
```

---

## 8. If something goes wrong

* **NaN at high energies**: the $\omega/2\omega$ coverage limit
  (Section 4). Check `stack.usable_omega(...)`; restrict the grid or pick a
  substrate record with wider coverage.
* **A method raises "requires Lz > 0" or "requires h_2D > 0"**: call
  `CHI.supercell_height_SI()` first, and construct the loader with `h_2D=` if
  you need effective quantities or $n$, $k$. The errors say exactly what is
  missing.
* **Numbers changed after a database update**: you selected records by index;
  re-select by `source=` and record `database_version()`.
* **Two runs give very different raw $\chi^{(2)}$**: compare `Lz` first; the raw
  values scale with vacuum padding. Compare `SHG_sheet`, not `order2`.
* **z component isn't zero**: stop as the run itself needs fixing before any
  analysis.

---

## 9. References

1. **Song et al.**, *Opt. Express* **31**, 19746 (2023) - structure-factor model implemented in `Stack`.
   [doi.org/10.1364/OE.486719](https://doi.org/10.1364/OE.486719)
2. **Cheng, Sipe, Vermeulen & Guo**, *J. Phys. Photonics* **1**, 015002 (2019) - the sheet-optics framework the structure factor builds on.
   [doi.org/10.1088/2515-7647/aaeadb](https://doi.org/10.1088/2515-7647/aaeadb)
3. **Woodward et al.**, *2D Mater.* **4**, 011006 (2017) - strict-SI single-interface formula implemented in `WoodwardModel`.
   [doi.org/10.1088/2053-1583/4/1/011006](https://doi.org/10.1088/2053-1583/4/1/011006)
4. **Clark et al.**, *Phys. Rev. B* **90**, 121409(R) (2014) - bulk thin-slab reference formula implemented in `ClarkModel`.
   [doi.org/10.1103/PhysRevB.90.121409](https://doi.org/10.1103/PhysRevB.90.121409)
5. **Bloembergen & Pershan**, *Phys. Rev.* **128**, 606 (1962) - the original theory of radiation from a nonlinear sheet.
   [doi.org/10.1103/PhysRev.128.606](https://doi.org/10.1103/PhysRev.128.606)
6. **Butcher & Cotter**, *The Elements of Nonlinear Optics*, Cambridge University Press (1990) - textbook form of the bulk formula (Eq. 7.27).
   [doi.org/10.1017/CBO9781139167994](https://doi.org/10.1017/CBO9781139167994)
