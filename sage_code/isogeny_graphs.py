"""isogeny_graphs.py

    Here we check for "unrecorded" isogenies; see Section 9 of the paper

    ====================================================================

    This file is part of Quadratic Kenku Solver.

    Copyright (C) 2022 Barinder S. Banwait, Filip Najman, and Oana
    Padurariu

    Quadratic Kenku Solver is free software: you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation, either version 3 of
    the License, or any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <https://www.gnu.org/licenses/>.

    ====================================================================

"""

from sage.all import EllipticCurve, Matrix, pari, gp
from sage.interfaces.quit import expect_quitall
from utils import GENUS_ONE_LIST, GENUS_ZERO_LIST
import logging

from timeout import timeout, TimeoutError

logger = logging.getLogger(__name__)
ISOGENY_CLASS_TIMEOUT_S = 30


def isogeny_degrees(j, K=None):
    """This function was suggested to us by John Cremona - thanks John!"""
    if K is None:
        K = j.parent()
    E = EllipticCurve(j=K(j))
    C = E.isogeny_class()
    return set(C.matrix()[C.index(E)])


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


@timeout(ISOGENY_CLASS_TIMEOUT_S)
def attempt_gp_comp(j, K, d):
    try:
        L,M = isogeny_class_via_gp(j, K, d)
        logger.info("PARI/GP computation worked!! :)")
        return L, M
    except TimeoutError:
        # now we really can't do anything more
        logger.info("PARI/GP took too long; assuming no unrecorded isogenies "
        "here, but you should check this directly in GP with the script in "
        "the `gp_code` folder."
        )
        return None, None
    except Exception:
        # Again, if PARI/GP fails, then nothing more that can be done here
        # except to ask the user to verify this directly in GP
        # because it's likely to be an issue with the Sage/GP interface 
        logger.info("PARI/GP failed; assuming no unrecorded isogenies "
        "here, but you should check this directly in GP with the script in "
        "the `gp_code` folder."
        )
        return None, None

@timeout(ISOGENY_CLASS_TIMEOUT_S)
def timed_isogeny_class(E):
    return E.isogeny_class()


def isogeny_class_via_sage(j, K, d):
    E = EllipticCurve(j=K(j))
    logger.debug(f"Constructing isogeny graph with j-invariant {j} ...")
    try:
        C = timed_isogeny_class(E)
        logger.debug("Done.")
        return [F.j_invariant() for F in C] , C.matrix()
    except Exception as err_msg:
        logger.warning(f"Isogeny graph computation failed with message: {err_msg}. "
        "We will now attempt the same computation in PARI/GP ... "
        )
        L,M = attempt_gp_comp(j, K, d)
        return L,M


def isogeny_class_via_pari(j, K):
    """This was my first attempt, but unfortunately doesn't work for all
    j-invariants for an unknown reason"""
    E = EllipticCurve(j=K(j))
    Epari = pari(E)
    L, M = Epari.ellisomat(1)
    jInvs = [EllipticCurve(K, list(eRep)).j_invariant() for eRep in L]
    return jInvs, Matrix(M)


def unrecorded_isogenies(K, my_js, d, z=None, cm=False):
    """This is a wrapper for the above function, to obtain the "unrecorded isogenies"
    (in the sense of Mazur) from the j-invariants identified above"""

    isog_classes_j_invs = []
    isog_mats = []
    failed_dict = {}
    for j in my_js:
        jInvs, M = isogeny_class_via_sage(j, K, d)
        if jInvs is not None:
            isog_classes_j_invs.append(jInvs)
            isog_mats.append(M)
        else:
            # means we got a failure
            failed_dict[j] = z if z else 0


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

    return unrecorded_isogeny_degrees, failed_dict
