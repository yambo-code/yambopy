"""This file contains some functions useful in in the different modules of qepy"""

def fortran_bool(boolean):
    return {True:'.true.',False:'.false.'}[boolean]
