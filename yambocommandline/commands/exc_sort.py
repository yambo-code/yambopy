#
# Authors: MN
#
import argparse
import os
from yambopy.dbs.excitondb import YamboExcitonDB
from yambopy.dbs.latticedb import YamboLatticeDB
#
def run_exc_sort(args):
    parser = argparse.ArgumentParser(description="Write the sorted energies and intensities to a file.")
    #
    parser.add_argument("--save", type=str, default=".", help="Path to the SAVE directory (location of ns.db1, as in PATH/SAVE). Default current directory.")
    parser.add_argument("-J", "--jobdir", type=str, default="SAVE", metavar="DIR", help="BSE job directory (as in PATH/DIR). Default SAVE")
    parser.add_argument("--iqpt", type=int, default=1, help="Q-point index. Default 1")
    args = parser.parse_args(args)

    nsfile = os.path.join(args.save, "ns.db1")
    lattice = YamboLatticeDB.from_db_file(nsfile)
    #
    filename = f"ndb.BS_diago_Q{args.iqpt}"
    excdb = YamboExcitonDB.from_db_file(lattice, filename=filename,
                                        folder=args.jobdir,
                                        Load_WF=False, neigs=-1)
    excdb.write_sorted(prefix=f"o-{args.jobdir.split('/')[-1]}.exc_qpt{args.iqpt}_sorted")
