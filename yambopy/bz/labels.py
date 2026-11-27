"""
labels.py
=========
Special-point labels of Brillouin-zone paths: parsing of path strings and
rendering of labels as plain text, LaTeX (matplotlib) or Unicode.

Grammar
-------
A path string is a comma-separated list of sections (a comma marks a
discontinuity); a section is a sequence of labels, optionally separated
by spaces or hyphens for readability::

    'GMKGALHA,LM,KH'   'G-M-K-G'   'Gamma M K Gamma'   "GKK'G"

A label is::

    base digits? "'"*

* ``base``: one capital Latin letter (``A``-``Z``), or a capitalised Greek
  name from :data:`GREEK` (``Gamma``, ``Sigma``, ...).  Since a Latin base
  is a single letter, a capital followed by lowercase letters can only be
  a Greek name, so the split is unambiguous.
* ``digits``: a subscript (``X1`` is X₁, ``Sigma1`` is Σ₁).
* ``'``: primes (``K'`` is K′).

``G`` is the name of Γ (as in ASE and the lattice data); ``Gamma`` is an
accepted alias and is stored as ``G``.  Greek capitals identical to Latin
ones (Alpha, Beta, Rho, ...) are not names: use the Latin letter.

Rendering (:func:`format_label`)
--------------------------------
=========  =========  ================  ==========  ================
fmt        ``G``      ``Sigma1``        ``K'``      merged ``A, L``
=========  =========  ================  ==========  ================
plain      ``G``      ``Sigma1``        ``K'``      ``A — L``
latex      ``$\\Gamma$``  ``$\\Sigma_{1}$``  ``$K^{\\prime}$``  ``$A|L$``
unicode    ``Γ``      ``Σ₁``            ``K′``      ``A|L``
=========  =========  ================  ==========  ================
"""

import re

__all__ = [
    "GREEK",
    "MERGE_SEPARATOR",
    "parse_path",
    "parse_label",
    "normalise_label",
    "label_to_latex",
    "label_to_unicode",
    "format_label",
]

#: Greek capitals usable as label bases (those distinct from Latin ones).
GREEK = {
    "Gamma": "Γ",
    "Delta": "Δ",
    "Theta": "Θ",
    "Lambda": "Λ",
    "Xi": "Ξ",
    "Pi": "Π",
    "Sigma": "Σ",
    "Upsilon": "Υ",
    "Phi": "Φ",
    "Psi": "Ψ",
    "Omega": "Ω",
}

#: Separator joining the end of a path section and the start of the next
#: in plain merged labels (e.g. ``'A — L'``).
MERGE_SEPARATOR = " — "

#: Plain names that are stored under another name.
_ALIASES = {"Gamma": "G"}

#: Latin letters rendered as Greek.
_LATIN_AS_GREEK = {"G": "Gamma"}

# Greek names are tried first (longest first, so that e.g. "Phi" is not
# read as "P"), then a single Latin capital.
_BASE = "|".join(sorted(GREEK, key=len, reverse=True)) + "|[A-Z]"
_LABEL_RE = re.compile(rf"(?P<base>{_BASE})(?P<sub>[0-9]*)(?P<primes>'*)")
_SEPARATORS = " \t-"

_SUBSCRIPTS = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

_LABEL_HELP = (
    "Labels are one capital letter or a capitalised Greek name "
    f"({', '.join(GREEK)}), optionally followed by digits and primes, "
    "e.g. G, X1, Sigma1, K'."
)


def parse_label(label: str) -> tuple[str, str, int]:
    """
    Split a label into ``(base, subscript, number of primes)``, e.g.
    ``"Sigma1"`` -> ``("Sigma", "1", 0)`` and ``"K'"`` -> ``("K", "", 1)``.

    Raises
    ------
    ValueError
        If *label* does not follow the label grammar.
    """
    match = _LABEL_RE.fullmatch(label)
    if match is None:
        raise ValueError(f"Invalid special-point label {label!r}. {_LABEL_HELP}")
    return match["base"], match["sub"], len(match["primes"])


def normalise_label(label: str) -> str:
    """Validate *label* and return its stored name (``Gamma...`` -> ``G...``)."""
    base, sub, primes = parse_label(label)
    return _ALIASES.get(base, base) + sub + "'" * primes


def parse_path(path: str) -> list[list[str]]:
    """
    Split a path string into sections of normalised labels, e.g.
    ``"GMKG,A-L"`` -> ``[['G', 'M', 'K', 'G'], ['A', 'L']]``.

    Empty sections (e.g. from a trailing comma) are skipped.

    Raises
    ------
    ValueError
        At the first character that does not start a valid label.
    """
    sections: list[list[str]] = []
    for section in path.split(","):
        labels: list[str] = []
        pos = 0
        while pos < len(section):
            if section[pos] in _SEPARATORS:
                pos += 1
                continue
            match = _LABEL_RE.match(section, pos)
            if match is None:
                raise ValueError(
                    f"Invalid path {path!r}: cannot read a label at {section[pos:]!r}. "
                    + _LABEL_HELP
                )
            labels.append(normalise_label(match.group()))
            pos = match.end()
        if labels:
            sections.append(labels)
    return sections


def label_to_latex(label: str, math: bool = True) -> str:
    """
    LaTeX for a label, e.g. ``"G"`` -> ``"$\\Gamma$"``, ``"X1"`` ->
    ``"$X_{1}$"``, ``"K'"`` -> ``"$K^{\\prime}$"``.

    With ``math=False`` the ``$`` delimiters are omitted.
    """
    base, sub, primes = parse_label(label)
    base = _LATIN_AS_GREEK.get(base, base)
    tex = "\\" + base if base in GREEK else base
    if sub:
        tex += "_{" + sub + "}"
    if primes:
        tex += "^{" + "\\prime" * primes + "}"
    return f"${tex}$" if math else tex


def label_to_unicode(label: str) -> str:
    """
    Unicode for a label, e.g. ``"G"`` -> ``"Γ"``, ``"Sigma1"`` -> ``"Σ₁"``,
    ``"K'"`` -> ``"K′"``.
    """
    base, sub, primes = parse_label(label)
    base = _LATIN_AS_GREEK.get(base, base)
    return GREEK.get(base, base) + sub.translate(_SUBSCRIPTS) + "′" * primes


def format_label(label: str, fmt: str = "plain") -> str:
    """
    Render a plain label, or a merged one (``'A — L'``, see
    :data:`MERGE_SEPARATOR`), in format *fmt*: ``'plain'`` (unchanged),
    ``'latex'`` (``'$A|L$'``) or ``'unicode'`` (``'A|L'``).

    Labels already written in LaTeX (containing ``$``) are returned
    unchanged.
    """
    if fmt == "plain" or "$" in label:
        return label
    parts = label.split(MERGE_SEPARATOR)
    if fmt == "latex":
        return "$" + "|".join(label_to_latex(p, math=False) for p in parts) + "$"
    if fmt == "unicode":
        return "|".join(label_to_unicode(p) for p in parts)
    raise ValueError(f"Unknown label format {fmt!r}: use 'plain', 'latex' or 'unicode'.")
