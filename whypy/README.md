##Development plan:

- Switch to a new structure in which objects are exposed (in the style of the yambo code), with separation of io, plot, data, etc from the individual classes. Objects could be scf/, nscf/, phonons/, bse/, screening/, gw/, etc. Distinction between QE and Yambo is hidden in io and data.
- First step: the plot directory, containing all plotting scripts. These will replace the plotting function within the various classes of yambopy/dbs and qepy, which will be removed in time.
