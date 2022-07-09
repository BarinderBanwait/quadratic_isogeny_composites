"""quadratic_kenku_solver.py

This generates entries for the table in the appendix

"""
import json
import logging
from sage.all import QuadraticField, cm_j_invariants, magma_free, EllipticCurve
from utils import minimally_finite_fast, GENUS_ONE_LIST, HYPERELLIPTIC_VALUES
from large_possible_isogeny_primes import LPIP
from isogeny_graphs import unrecorded_isogenies
from hyperelliptic_verifs import try_trbovic_filter
from non_hyperelliptic_verifs import (
    try_oezman_sieve,
    check_mwgp_same_minus,
)
from utils import check_mwgp_same_plus
from functools import reduce

QUADRATIC_POINTS_DATA_PATH = "quadratic_points_catalogue.json"

with open(QUADRATIC_POINTS_DATA_PATH, "r") as qdpts_dat_file:
    qdpts_dat = json.load(qdpts_dat_file)


def reducer(accumulator, element):
    for key, value in element.items():
        accumulator[key] = accumulator.get(key, 0) + value
    return accumulator


def format_elliptic_count_magma_function(d):

    nl = chr(10)

    with open("EllipticCount_magma_function.txt", "r") as magma_file:
        the_lines = magma_file.read().splitlines()

    return nl.join(the_lines) + nl + f"EllipticCount({d});"


def format_preimages_magma_function(d):

    nl = chr(10)

    with open("ComputePreimages_magma_function.txt", "r") as magma_file:
        the_lines = magma_file.read().splitlines()

    return nl.join(the_lines) + nl + f"main({d});"



def unique_j_inv_count(j_invs_str, K_gen):
    """This function takes input a list of strings representing
    j-invariants, and returns the number of unique ones. It is slightly
    more complicated than just taking len, since the representation is
    only up to Galois conjugacy, and it is in principle possible that
    some of them are Q-rational
    """
    if not j_invs_str:
        return 0, []
    j_invs_unpacked = []
    for j_inv_str_rep in j_invs_str:
        a, b = eval(j_inv_str_rep)
        if b == 0:
            j_invs_unpacked.append(a)
        else:
            j1 = a + b * K_gen
            j2 = a - b * K_gen
            j_invs_unpacked.append(j1)
            j_invs_unpacked.append(j2)
    return len(j_invs_unpacked), j_invs_unpacked


def process_hyperelliptic(d, K_gen, hyperelliptic_vals):

    output_dict = {}
    K = K_gen.parent()

    for z in hyperelliptic_vals:
        data_this_z = qdpts_dat[str(z)]
        if z != 37:
            if try_trbovic_filter(d, z):
                if try_oezman_sieve(d, z):
                    raise NotImplementedError(f"d = {d}, z={z}")
            j_invs_str = data_this_z["non_cm_points"].get(str(d), [])
            isog_count, j_inv_list = unique_j_inv_count(j_invs_str, K_gen)
            logging.debug(f"Starting isogeny class computation for j-invariants on X0({z})")
            unrecorded_isogenies_dict = unrecorded_isogenies(K, j_inv_list, d, z=z)
            logging.debug(f"DONE isogeny class computation for j-invariants on X0({z})!")
            # The following will only update the dictionary with unrecorded isogenies,
            # that is, with degrees which are multiples of z. Updating z itself
            # is done right at the end of the for loop
            output_dict = {**output_dict, **unrecorded_isogenies_dict}
        else:
            # 37 requires special handling
            if EllipticCurve("37b1").quadratic_twist(d).analytic_rank() == 0:
                # import pdb; pdb.set_trace()
                if eval(str(magma_free.eval(format_preimages_magma_function(d)))) == 4:
                    # then we only have the two Q-rational j-invariants
                    isog_count = 2
                    j_inv_list = [-162677523113838677, -9317]
                    logging.debug(f"Starting isogeny class computation for j-invariants on X0({z})")
                    unrecorded_isogenies_dict = unrecorded_isogenies(K, j_inv_list, d, z=z)
                    logging.debug(f"DONE isogeny class computation for j-invariants on X0({z})!")
                    output_dict = {**output_dict, **unrecorded_isogenies_dict}
                else:
                    raise NotImplementedError
            elif check_mwgp_same_plus(37, d):
                if not try_oezman_sieve(d, 37):
                    isog_count = 2
                    j_inv_list = [-162677523113838677, -9317]
                    logging.debug(f"Starting isogeny class computation for j-invariants on X0({z})")
                    unrecorded_isogenies_dict = unrecorded_isogenies(K, j_inv_list, d, z=z)
                    logging.debug(f"DONE isogeny class computation for j-invariants on X0({z})!")
                    output_dict = {**output_dict, **unrecorded_isogenies_dict}
                else:
                    raise NotImplementedError
            else:
                raise NotImplementedError

        output_dict[z] = isog_count

    return output_dict


def process_non_hyperelliptic(d, K_gen, non_hyperelliptic_vals):

    non_hyperelliptic = non_hyperelliptic_vals.copy()
    non_hyperelliptic.remove(
        163
    )  # we know they are all CM, so will be dealt with hereafter

    output_dict = {}
    K = K_gen.parent()

    for z in non_hyperelliptic:
        if str(z) in qdpts_dat:
            # means that the exceptional points have been determined
            data_this_z = qdpts_dat[str(z)]
            if not data_this_z["is_complete"]:
                # means that there are infinitely many quadratic
                # points. This would clearly be a problem; however we
                # try to rule out the infinitely many non-exceptional points
                # via the Ozman sieve and method of appendix.
                if try_oezman_sieve(d, z):
                    # This means that the Ozman sieve failed to rule out
                    # non-exceptional points. In this case we try the
                    # method of the appendix
                    if not check_mwgp_same_minus(z, d):
                        # If we ever get here, that means that BOTH Ozman sieve
                        # and method of appendix failed, and we must stop
                        raise NotImplementedError(
                            "The Ozman sieve failed. There are other things"
                            "which can be tried (like quadratic Chabauty) which might be added eventually."
                        )

            # if we get to this point, then the quadratic points catalogue
            # contains all the information we need. Thus we get the non CM
            # j-invariants, and construct the isogeny graphs to get the
            # unrecorded isogenies
            j_invs_str = data_this_z["non_cm_points"].get(str(d), [])
            isog_count, j_inv_list = unique_j_inv_count(j_invs_str, K_gen)
            logging.debug(f"Starting isogeny class computation for j-invariants on X0({z})")
            unrecorded_isogenies_dict = unrecorded_isogenies(K, j_inv_list, d, z=z)
            logging.debug(f"DONE isogeny class computation for j-invariants on X0({z})!")
            # The following will only update the dictionary with unrecorded isogenies,
            # that is, with degrees which are multiples of z. Updating z itself
            # is done right at the end of the for loop
            output_dict = {**output_dict, **unrecorded_isogenies_dict}
        else:
            # this means that the exceptional points have not been determined. In this
            # case we can still try the method of the Appendix

            if check_mwgp_same_minus(z, d):
                # this means that the method of the appendix worked, i.e., that either the
                # Q-sqrt-d points are all rational, or CM. Since CM points will be dealt with
                # hereafter, we have that they are all rational. Since z here is not in the
                # catalogue, and is a non-hyperelliptic value, we have that z must be 163.
                # But 163 was removed from the list above; therefore there are no rational
                # points on X0(z), so no isogenies
                isog_count = 0
            elif check_mwgp_same_plus(z, d):
                # If the method of the appendix failed, then we can try the "no growth in plus part"
                # method. If this succeeds, then it tells us that any Q-sqrt-d point must be
                # non-exceptional, and therefore we can try the Ozman sieve to rule it out
                if not try_oezman_sieve(d, z):
                    # Ozman sieve worked, so again we have no isogenies
                    isog_count = 0
                else:
                    raise NotImplementedError
            else:
                raise NotImplementedError

        output_dict[z] = isog_count

    return output_dict


def quadratic_kenku_solver(d):
    """Attempts to determine all cyclic isogenies of elliptic curves
    over Q-sqrt-d. Returns string of data for the table in the appendix
    if successful, returns None otherwise"""

    # First get the list of Minimally Finite values

    if not d in LPIP:
        raise ValueError(
            f"no large possible isogeny prime data at {d}. Please add to that file first."
        )

    K = QuadraticField(d, 'K_gen')
    K_gen = K.gen()
    minimally_finite = minimally_finite_fast(d)
    minimally_finite.extend(LPIP[d])
    minimally_finite.sort()
    logging.info(f"minimally finite values are {minimally_finite}")

    infinite_col = [N for N in GENUS_ONE_LIST if not N in minimally_finite]
    logging.info(f"infinite_col = {infinite_col}")

    my_elliptic_vals = []
    my_hyperelliptic_vals = []
    my_non_hyperelliptic_vals = []

    for mf_val in minimally_finite:
        if mf_val in GENUS_ONE_LIST:
            my_elliptic_vals.append(mf_val)
        elif mf_val in HYPERELLIPTIC_VALUES:
            my_hyperelliptic_vals.append(mf_val)
        else:
            my_non_hyperelliptic_vals.append(mf_val)

    logging.info(f"my_elliptic_vals = {my_elliptic_vals}")
    logging.info(f"my_hyperelliptic_vals = {my_hyperelliptic_vals}")
    logging.info(f"my_non_hyperelliptic_vals = {my_non_hyperelliptic_vals}")

    logging.info("Starting elliptic values ...")
    elliptic_jInv_magma_str = str(magma_free.eval(format_elliptic_count_magma_function(d)))
    elliptic_jInv_magma_str = elliptic_jInv_magma_str.replace("<", "(")
    elliptic_jInv_magma_str = elliptic_jInv_magma_str.replace(">", ")")
    elliptic_jInv_dict = dict(eval(elliptic_jInv_magma_str))
    elliptic_jInv_dict = {k:[K(a) for a in elliptic_jInv_dict[k] ] for k in elliptic_jInv_dict}

    elliptic_count_dict = {k:len(elliptic_jInv_dict[k]) for k in elliptic_jInv_dict}

    elliptic_unrecorded_isogenies_dict = {}
    for k in elliptic_count_dict:
        unrecorded_isogenies_dict = unrecorded_isogenies(K, elliptic_jInv_dict[k], d, z=k)
        elliptic_unrecorded_isogenies_dict = {**elliptic_unrecorded_isogenies_dict, **unrecorded_isogenies_dict}

    elliptic_count_dict = {**elliptic_count_dict, **elliptic_unrecorded_isogenies_dict}
    logging.info("Done elliptic values!")

    logging.info("Starting hyperelliptic values ...")
    hyperelliptic_count_dict = process_hyperelliptic(d, K_gen, my_hyperelliptic_vals)
    logging.info("Done hyperelliptic values!")

    logging.info("Starting non-hyperelliptic values ...")
    non_hyperelliptic_count_dict = process_non_hyperelliptic(
        d, K_gen, my_non_hyperelliptic_vals
    )
    logging.info("Done non-hyperelliptic values!")

    print(f"hyperelliptic_count = {hyperelliptic_count_dict}")
    print(f"non_hyperelliptic_count = {non_hyperelliptic_count_dict}")

    cm_jinvs = cm_j_invariants(K)

    logging.debug("Starting CM computation ...")
    cm_isogenies_dict = unrecorded_isogenies(K, cm_jinvs, d, cm=True)
    logging.debug("Computed CM isogeny classes!")
    all_dicts = [
        elliptic_count_dict,
        hyperelliptic_count_dict,
        non_hyperelliptic_count_dict,
        cm_isogenies_dict,
    ]
    ans = reduce(reducer, all_dicts, {})
    ans = {k: ans[k] for k in ans if ans[k] != 0}
    return ans


if __name__ == "__main__":
    verbose = True
    loglevel = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s",
        datefmt="%H:%M:%S",
        level=loglevel,
    )
    logging.debug("Debugging level for log messages set.")
    ans = quadratic_kenku_solver(213)
    logging.info(f"final answer is {ans}")
