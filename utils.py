"""utils.py

Some helpful functions
"""

def split_cartan_genus(p):
    """Computes the genus of X_s(p) from Imin Chen's "The Jacobians of
    non-split Cartan modular curves"""

    return (1/24) * (p^2 - 8*p + 11 - 4 * kronecker_symbol(-3,p))
