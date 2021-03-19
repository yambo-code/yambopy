##Development plan:

- Switch to a new structure in which objects are exposed (in the style of the yambo code), with separation of io, plot, data, job, etc from the individual classes. 

- "Job" objects could be scf, nscf, phonons, bse, screening, gw, etc. Distinction between QE and Yambo is hidden in io and data.

- "Data" objects could be kpoints, band structure, BZ data, spectral function, excitons, electrons, QPs, etc.

- Including database structure to keep track and organise calculations automatically

- First step: the plot directory, containing all plotting scripts and a plotting driver. These will replace the plotting function within the various classes of yambopy/dbs and qepy, which will be removed in time.
