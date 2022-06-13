"""isogeny_graphs.py

    Here we check for "unrecorded" isogenies; see Section 9 of the paper

"""

# We collect the j-invariants whose isogeny graphs we need to construct

my_js_213 = [-1159088625/2097152, -189613868625/128,3375/2,-140625/8, -884736, -882216989/131072, -297756989/2,-3375, 16581375, -32768, -24729001, -121, -9317, -162677523113838677, -884736000, -147197952000, -262537412640768000]
my_js_438 = [-884736, -882216989/131072, -297756989/2, -32768, -24729001, -121, -9317, -162677523113838677, -884736000, -147197952000, -262537412640768000]

# Obtainining the unrecorded isogeny degrees is then done as follows

def isogeny_degrees(j, K=None):
    """This function was suggested to us by John Cremona - thanks John!"""
    if K is None:
        K = j.parent()
    E = EllipticCurve(j=K(j))
    C = E.isogeny_class()
    return set(C.matrix()[C.index(E)])


def get_unrecorded_isogenies(d, my_js):
    """This is a wrapper for the above function, to obtain the "unrecorded isogenies"
    (in the sense of Mazur) from the j-invariants identified above"""
    K.<a> = QuadraticField(d)
    unrecorded_isogeny_degrees = set()
    for j in my_js:
        unrecorded_isogeny_degrees = unrecorded_isogeny_degrees.union(isogeny_degrees(K(j)))

    return unrecorded_isogeny_degrees

# Then run the following to verify the main claim of Section 9 of the paper

# get_unrecorded_isogenies(213, my_js_213)
# get_unrecorded_isogenies(438, my_js_438)