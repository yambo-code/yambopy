import numpy as np 

def find_degeneracy_evs(eigenvalues, atol=1e-3, rtol=1e-3):
    """
    Identify sets of degenerate eigenstates based on eigenvalues.

    Parameters
    ----------
    eigenvalues : array-like
        A sorted list or array of eigenvalues in increasing order.
    atol : float, optional
        Absolute tolerance for degeneracy comparison. Default is 1e-3.
    rtol : float, optional
        Relative tolerance for degeneracy comparison. Default is 1e-3.

    Returns
    -------
    list of lists
        A list where each sublist contains the indices of degenerate eigenstates.

    Raises
    ------
    ValueError
        If `eigenvalues` is empty or not a valid array.
    """
    # Input validation
    eigenvalues = np.asarray(eigenvalues)
    if eigenvalues.size == 0:
        raise ValueError("Input eigenvalues must not be empty.")
    if atol < 0 or rtol < 0:
        raise ValueError("Tolerances `atol` and `rtol` must be non-negative.")

    # Compute differences between consecutive eigenvalues
    diffs = np.diff(eigenvalues)
    tolerance = atol + rtol * np.abs(eigenvalues[:-1])

    # Identify where the differences exceed the tolerance
    split_indices = np.where(diffs > tolerance)[0]

    # Group indices of degenerate states
    degen_sets = np.split(np.arange(len(eigenvalues)), split_indices + 1)

    # Convert numpy arrays to lists for consistency
    degen_sets = [group for group in degen_sets]

    return degen_sets


if __name__ == '__main__':
    eigenvalues = [1.0, 1.001, 1.002, 2.0, 2.001, 3.0]
    degen_sets = find_degeneracy_evs(eigenvalues, atol=1e-3, rtol=1e-3)
    print(degen_sets)
