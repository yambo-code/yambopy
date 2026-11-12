#
# Authors: MN
#
from yambopy.bse.exciton_irreps import compute_exc_rep
import argparse
#
def run_exc_irrep(args):
    """
    Compute the irreducible representation (irrep) of excitonic states.

    This is a command-line interface helper that parses input arguments
    and invokes :func:`compute_exc_rep` to analyze exciton symmetry
    properties from a BSE (Bethe–Salpeter Equation) calculation.

    Parameters
    ----------
    args : list[str]
        Command-line arguments list (e.g. ``sys.argv[1:]``). The supported options are:

        * ``--path`` : str, optional
          Path to the main calculation directory (default: ``"."``).

        * ``-J`` / ``--jobdir`` : str, optional
          Subdirectory containing the BSE results (default: ``"SAVE"``).

        * ``--iqpt`` : int, optional
          Index of the transferred momentum (q-point) to analyze (default: ``1``).

        * ``--nstates`` : int, optional
          Number of exciton states to process (default: ``1``).

        * ``--degen_tol`` : float, optional
          Energy tolerance used to detect degenerate excitons, in eV (default: ``1e-2``).
          Total tolerence = degen_rtol*Energy_in_eV + degen_tol

        * ``--degen_rtol`` : float, optional
          Relative tolerance used to detect degenerate excitons (default: ``1e-3``).
          Total tolerence = degen_rtol*Energy_in_eV + degen_tol

        * ``--sym_tol`` : float, optional
          Numerical tolerance for evaluating symmetry operations (default: ``1e-2``).

        * ``--use_save_symm`` : bool, optional
          If True, use save symmetries instead of spglib ones. For non-collinear magnets, it is automatically
          set to True,

    Notes
    -----
    This function only handles argument parsing and dispatching.
    """
    parser = argparse.ArgumentParser(description="Compute exciton representation of excitonic states.")
    #
    parser.add_argument("--path", type=str, default=".", help="Path to the calculation directory. Default current directory.")
    parser.add_argument("-J", "--jobdir", type=str, default="SAVE", metavar="DIR", help="BSE job directory. Default SAVE")
    parser.add_argument("--iqpt", type=int, default=1, help="Q-point index. Default 1")
    parser.add_argument("--nstates", type=int, default=1, help="Number of exciton states. Default 1")
    parser.add_argument("--degen_tol", type=float, default=1e-2, help="Absolute tolerance for degeneracy. " +
                        "Default 0.01 eV. Total_tolerence = degen_rtol*Energy_in_eV + degen_tol.")
    parser.add_argument("--degen_rtol", type=float, default=1e-3, help="Relative tolerance for degeneracy. "+
                        "Default 0.001. Total_tolerence = degen_rtol*Energy_in_eV + degen_tol.")
    parser.add_argument("--sym_tol", type=float, default=1e-2, help="Tolerance for Symmetry operations. Default 0.01")
    parser.add_argument("--use_save_symm", action="store_true", help="use save symmetries instead of spglib ones. "
                        "For non-collinear magnets, it is automatically set to True")
    #
    args = parser.parse_args(args)
    compute_exc_rep(path=args.path,bse_dir=args.jobdir, iqpt=args.iqpt,
                    nstates=args.nstates, degen_tol=args.degen_tol,degen_rtol=args.degen_rtol,
                    symm_tol=args.sym_tol,use_save_symmetries=args.use_save_symm)
