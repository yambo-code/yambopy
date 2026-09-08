def latex(s: str) -> str:
    """
    Print to terminal Greek letters
    written in latex math mode
    """
    import re
    LATEX_GREEK = {
        r'\alpha': 'α', \
        r'\beta': 'β', \
        r'\gamma': 'γ', \
        r'\delta': 'δ',\
        r'\epsilon': 'ε',\
        r'\zeta': 'ζ',\
        r'\eta': 'η',\
        r'\theta': 'θ',\
        r'\iota': 'ι',\
        r'\kappa': 'κ',\
        r'\lambda': 'λ',\
        r'\mu': 'μ',\
        r'\nu': 'ν',\
        r'\xi': 'ξ',\
        r'\omicron': 'ο',\
        r'\pi': 'π',\
        r'\rho': 'ρ',\
        r'\sigma': 'σ',\
        r'\tau': 'τ',\
        r'\upsilon': 'υ',\
        r'\phi': 'φ',\
        r'\chi': 'χ',\
        r'\psi': 'ψ',\
        r'\omega': 'ω',
        r'\Gamma': 'Γ',\
        r'\Delta': 'Δ',\
        r'\Theta': 'Θ',\
        r'\Lambda': 'Λ',\
        r'\Xi': 'Ξ',\
        r'\Pi': 'Π',\
        r'\Sigma': 'Σ',\
        r'\Upsilon': 'Υ',\
        r'\Phi': 'Φ',\
        r'\Psi': 'Ψ',\
        r'\Omega': 'Ω',
    }

    # remove dollars
    s = re.sub(r'\$(.*?)\$', r'\1', s)
    # replace latex commands
    for symbol, unicode in LATEX_GREEK.items():
        s = s.replace(symbol, unicode)

    return s
