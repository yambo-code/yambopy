import numpy as np

def find_degeneracy_evs(eigenvalues, atol=1e-3, rtol=1e-3):
    """
    Identify sets of degenerate eigenstates based on eigenvalues.

    Parameters
    ----------
    eigenvalues : array-like
        A list or array of eigenvalues
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
    eigenvalues = np.asarray(eigenvalues)
    
    if eigenvalues.size == 0:
        raise ValueError("Input eigenvalues must not be empty.")
    if atol < 0 or rtol < 0:
        raise ValueError("Tolerances `atol` and `rtol` must be non-negative.")
    
    # Sort eigenvalues and get sorted indices
    idx_sorted = np.argsort(eigenvalues)
    eigenvalues = eigenvalues[idx_sorted]
    
    # Compute differences between consecutive eigenvalues
    diffs = np.diff(eigenvalues)
    tolerance = atol + rtol * np.abs(eigenvalues[:-1])
    
    # Identify where the differences exceed the tolerance
    split_indices = np.where(diffs > tolerance)[0]
    
    # Group indices of degenerate states
    degen_sets = np.split(idx_sorted, split_indices + 1)
    
    # Further split groups based on the mean of the group
    # NM : This is to ensure that the list is a Arithmetic progression
    # with d < tol
    final_degen_sets = []
    for group in degen_sets:
        if len(group) == 0:
            continue
        group_eigenvalues = eigenvalues[group]
        current_group = [group[0]]
        current_mean = group_eigenvalues[0]
        
        for i in range(1, len(group)):
            diff = np.abs(group_eigenvalues[i] - current_mean)
            tolerance = atol + rtol * np.abs(current_mean)
            
            if diff <= tolerance:
                current_group.append(group[i])
                # Update the mean of the current group incrementally
                current_mean = (current_mean * (len(current_group) - 1) + group_eigenvalues[i]) / len(current_group)
            else:
                final_degen_sets.append(current_group)
                current_group = [group[i]]
                current_mean = group_eigenvalues[i]
        
        # Append the last group
        if current_group:
            final_degen_sets.append(current_group)
    
    return final_degen_sets


if __name__ == '__main__':
    eigenvalues = [3.0001, 2.99976, 1.99999, 1.0, 1.001, 1.002, 2.0, 2.001, 3.0]
    degen_sets = find_degeneracy_evs(eigenvalues, atol=1e-3, rtol=1e-3)
    print(degen_sets)
