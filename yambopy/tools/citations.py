#
# Authors: MN
#
import atexit

# Global citation registry
_CITATIONS_USED = set()


def citation(ref: str):
    """
    Decorator that records a citation reference when the function is executed.

    This decorator may be applied multiple times to the same function,
    enabling stacked usage such as::

        @citation("A")
        @citation("B")
        def func():
            ...

    Each reference is recorded when the wrapped function is run. Duplicates
    across multiple calls are ignored.

    Parameters
    ----------
    ref : str
        A reference string identifying the source to cite, such as
        a publication title, DOI, or URL.

    Returns
    -------
    callable
        A decorated function that appends the citation to the global registry
        whenever invoked.
    """
    def wrapper(func):
        def inner(*args, **kwargs):
            _CITATIONS_USED.add(ref)
            return func(*args, **kwargs)

        # Preserve function identity for debugger, help(), and Sphinx autodoc
        inner.__name__ = func.__name__
        inner.__doc__ = func.__doc__
        inner.__qualname__ = func.__qualname__
        return inner
    return wrapper


@atexit.register
def print_citations():
    """
    Print a numbered citation list on interpreter shutdown.

    If any decorated functions were executed during the program run,
    this function prints a formatted list of unique citation references
    in sorted order.

    Notes
    -----
    This function is automatically registered via :mod:`atexit`.
    Users should not call this function directly.

    Examples
    --------
    Typical output:: Run this script.
    """
    if not _CITATIONS_USED:
        return

    print("\n=========================================================")
    print("Please cite the following references in case you wish")
    print("to acknowledge the work done by the Authors.")
    # Always print yambopy citation in the start.
    print("=========================================================")
    print("1) F. Paleari, et al. Yambopy. Zenodo (2025). doi: 10.5281/zenodo.15012963")
    for i, ref in enumerate(sorted(_CITATIONS_USED), start=1):
        print(f"{i+1}) {ref}")
    print("=========================================================")


#### test
if __name__ == "__main__":
    @citation("Ref A: Tenet by C. Nolan")
    def Tenet():
        return 1 + 1

    @citation("Ref B: Intersteller by C. Nolan")
    @citation("Ref C: Oppenheimer by C. Nolan")
    def Intersteller():
        return 1 + 1
    # call functcion
    for i in range (3) :
        Intersteller()
        Tenet()
