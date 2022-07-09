"""isogeny_graphs.py

    Here we check for "unrecorded" isogenies; see Section 9 of the paper

"""

from sage.all import EllipticCurve, Matrix, pari, gp
from sage.interfaces.quit import expect_quitall
from utils import GENUS_ONE_LIST, GENUS_ZERO_LIST
import logging

logger = logging.getLogger(__name__)

def isogeny_degrees(j, K=None):
    """This function was suggested to us by John Cremona - thanks John!"""
    if K is None:
        K = j.parent()
    E = EllipticCurve(j=K(j))
    C = E.isogeny_class()
    return set(C.matrix()[C.index(E)])


def isogeny_class_via_sage(j, K, d):

    E = EllipticCurve(j=K(j))
    C = E.isogeny_class()
    return [F.j_invariant() for F in C] , C.matrix()

def isogeny_class_via_pari(j, K):
    """This was my first attempt, but unfortunately doesn't work for all
    j-invariants for an unknown reason"""
    E = EllipticCurve(j=K(j))
    Epari = pari(E)
    L, M = Epari.ellisomat(1)
    jInvs = [EllipticCurve(K, list(eRep)).j_invariant() for eRep in L]
    return jInvs, Matrix(M)


def isogeny_class_via_gp(j, K, d):
    """The above function using the PARI C library functions worked for
    some j-invariants, but not for others. I have no idea why?!? Anyway,
    this function pacakges everything for computation directly in GP.
    """
    gp("default(parisize,256000000)");
    gp(f"K = nfinit(a^2 - {d})");
    j_as_list = list(K(j))
    assert len(j_as_list) == 2
    gp(f"myJ = Mod({j_as_list[0]} + {j_as_list[1]} * a, a ^ 2 - {d})");
    gp("v = ellfromj(myJ)");
    gp("E = ellinit(v, K)");
    logger.debug(f"Constructing isogeny graph with j-invariant {j} ...")
    L,M = gp("ellisomat(E, 1)");
    logger.debug("Done.")
    jInvs = [EllipticCurve(K, list(eRep)).j_invariant() for eRep in pari(L)]
    M_matrix = Matrix(pari(M))
    expect_quitall()
    return jInvs, M_matrix


def unrecorded_isogenies(K, my_js, d, z=None, cm=False):
    """This is a wrapper for the above function, to obtain the "unrecorded isogenies"
    (in the sense of Mazur) from the j-invariants identified above"""

    isog_classes_j_invs = []
    isog_mats = []

    for j in my_js:
        jInvs, M = isogeny_class_via_sage(j, K, d)
        isog_classes_j_invs.append(jInvs)
        isog_mats.append(M)

    if cm:
        desired_degrees = {
            a
            for M in isog_mats
            for v in M
            for a in v
            if a not in set(GENUS_ZERO_LIST).union(set(GENUS_ONE_LIST))
        }
    else:
        assert z is not None
        desired_degrees = {
            a for M in isog_mats for v in M for a in v if a % z == 0 and a > z
        }

    unrecorded_isogeny_degrees = {}
    for k in desired_degrees:
        j_invs_admitting_this_deg = set()
        for C, M in zip(isog_classes_j_invs, isog_mats):
            for i, v in enumerate(M):
                if k in v:
                    j_invs_admitting_this_deg.add(C[i])
        unrecorded_isogeny_degrees[k] = len(j_invs_admitting_this_deg)

    return unrecorded_isogeny_degrees


# Then run the following to verify the main claim of Section 9 of the paper

# unrecorded_isogenies(213, my_js_213)
# unrecorded_isogenies(438, my_js_438)

# K = nfinit(Pol(Vecrev([-1, -1, 1])))
# E = ellinit(
#     [
#         Pol(Vecrev([1, 0])),
#         Pol(Vecrev([1, 1])),
#         Pol(Vecrev([0, 1])),
#         Pol(Vecrev([0, 1])),
#         Pol(Vecrev([0, 0])),
#     ],
#     K,
# )
# [L, M] = ellisomat(E, 1)

# K = nfinit(a ^ 2 - 213)
# myJ = Mod(123 + 456 * a, a ^ 2 - 5)
# v = ellfromj(myJ)
# E = ellinit(v, K)
# ellisomat(E, 1)
# test = ellinit(L[1]).j
