#!/usr/bin/env python3
#
# Authors: MN
import argparse
import numpy as np
import os

from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.latticedb import YamboLatticeDB
from yambopy.dbs.wfdb import YamboWFDB


def run_plot_exc_wf_real_space(args):
    parser = argparse.ArgumentParser(description="Generate real-space " +
                                     "exciton wavefunction in Gaussian .cube file.")

    parser.add_argument("--path", type=str, default=".",
                        help="Calculation directory [PATH/SAVE] (default: Current working directory '.')")
    parser.add_argument("-J","--jobdir", type=str, default="SAVE", metavar="DIR",
                        help="BSE JOB directory [PATH/DIR] (default: SAVE)")
    parser.add_argument("--iqpt", type=int, default=1,help="Q-point index (default: 1)")
    parser.add_argument("--iexe", type=int, required=True,help="Exciton index to plot. (counts from 1)")
    # --iexe and --iqpt are not python indexing. i,e 1st item starts from 1.
    parser.add_argument("--supercell", nargs=3, type=int, default=[1, 1, 1],
                        help="Supercell dimensions (default: 1 1 1)")
    parser.add_argument("--wfc_cutoff", type=float, default=-1.0,
                        help="Wavefunction cutoff in Ry (default: −1 → full cutoff)")
    parser.add_argument("--degen_tol", type=float, default=0.01,
                        help="Degeneracy threshold in eV (default: 0.01)")
    # Mutually-exclusive: fix hole OR electron position
    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument("--hole", nargs=3, type=float,
                     help="Fix hole position in reduced units.")
    grp.add_argument("--electron", nargs=3, type=float,
                     help="Fix electron position in reduced units.")
    #
    parser.add_argument("--phase", action="store_true",
                        help="Enable phase plot (default: False). "
                        "Works only for non-degenerate excitons "
                        "and nspin=nspinor=1.")
    #
    parser.add_argument("--block_size", type=int, default=256,
                        help="Block size for computation (default: 256). "
                        "Lower only for very tight memory.")
    #
    parser.add_argument("--neig_load", type=int, default=100,
                        help="Number of extra eigevectors to load in memory (default: 100). "
                        "Only meaningful for TDA (non TDA, it is ignored and all are read). "
                        "Default 100 is fine and in most cases can be reduced, but if you are in deep "
                        "continum, we might have more accidental degenerates, in that case "
                        "we need to load more eigenvectors.")
    #
    args = parser.parse_args(args)
    #
    calc_path = args.path
    BSE_dir = args.jobdir
    iqpt = args.iqpt
    #
    if args.hole:
        fix_particle = "h"
        fixed_position = args.hole
    else:
        fix_particle = "e"
        fixed_position = args.electron
    #
    nsfile = os.path.join(calc_path, "SAVE", "ns.db1")
    #
    lattice = YamboLatticeDB.from_db_file(nsfile)
    #
    filename = f"ndb.BS_diago_Q{iqpt}"
    # NM : In the case of TDA, we donot want to load all the eigenvectors. 
    excdb = YamboExcitonDB.from_db_file(lattice, filename=filename,
                                        folder=BSE_dir,
                                        neigs=args.iexe + args.neig_load)
    #
    wfdb = YamboWFDB(path=calc_path,latdb=lattice,
                     bands_range=[np.min(excdb.table[:, 1]) - 1,np.max(excdb.table[:, 2])])
    #
    excdb.real_wf_to_cube(iexe=args.iexe-1,wfdb=wfdb,fixed_postion=fixed_position,
                          supercell=args.supercell,degen_tol=args.degen_tol,
                          wfcCutoffRy=args.wfc_cutoff,fix_particle=fix_particle,
                          phase=args.phase, block_size=args.block_size)



