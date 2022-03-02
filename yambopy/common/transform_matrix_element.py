# This file is part of yambopy
# Author: FP

class ExpandMatrixElement():
    """
    This is a generic symmetry expansion script.

    Input: 
      - System symmetries from YamboLatticeDB
      - Array of matrix elements on k/q-grids to be expanded from IBZ to BZ
      - Instruction to expand over k or over q

    Output:
      - Array of expanded matrix elements in the specified grid

    Example:

      :: O_nmR{k}^S{q} =  <n R{k}|O|mR{k}-S{q}>

      :: Expansion over q is

      :: O_nmR{k}^S{q} =  O_nm S^-1R{k}^q     (no TR)
      :: O_nmR{k}^S{q} = [O_nm S^-1R{k}^q]^*  (TR)

      :: Expansion over k is

      :: O_nmR{k}^S{q} =  O_nm k^R^-1S{q}     (no TR)
      :: O_nmR{k}^S{q} = [O_nm k^R^-1S{q}]^*  (TR)

    """
    #def __init__(mats_ibz,syms,space='q',TR=False):
        