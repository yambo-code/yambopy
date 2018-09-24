"""Duck-typing tests (taken from abipy)"""

def isstring(s):
    """True if s behaves like a string (duck typing test)."""
    try:
        s + " "
        return True
    except TypeError:
        return False

def isiter(s):
    try:
        iterator = iter(s)
        return True
    except TypeError:
        return False
