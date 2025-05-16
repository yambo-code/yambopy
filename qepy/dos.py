import numpy as np

def density_of_states(energies, energy_grid, width=0.1, broadening='gaussian'):
    """
    Calculate the Density of States (DOS) using either Gaussian or Lorentzian broadening.
    
    Parameters:
        energies (array-like): Array of energy eigenvalues.
        energy_grid (array-like): Grid of energy values where DOS will be evaluated.
        width (float): Broadening width (standard deviation for Gaussian, HWHM for Lorentzian).
        broadening (str): Type of broadening to use; options are 'gaussian' or 'lorentzian'.
    
    Returns:
        dos (np.ndarray): Density of States evaluated on the energy grid.
    """
    energies = np.array(energies)
    energy_grid = np.array(energy_grid)
    dos = np.zeros_like(energy_grid)

    # Choose broadening function
    if broadening.lower() == 'gaussian':
        # Gaussian broadening formula
        for e in energies:
            dos += np.exp(-((energy_grid - e)**2) / (2 * width**2)) / (width * np.sqrt(2 * np.pi))
    elif broadening.lower() == 'lorentzian':
        # Lorentzian broadening formula
        for e in energies:
            dos += (width / np.pi) / ((energy_grid - e)**2 + width**2)
    else:
        raise ValueError("Broadening must be 'gaussian' or 'lorentzian'")
    
    return dos
