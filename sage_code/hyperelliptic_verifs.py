"""hyperelliptic_verifs.py

This verifies the claims in Section 6 of the paper

"""

# We first show that, for all N as in the statement of Proposition 6.1,
# apart from N = 37, the Najman-Trbovic filter method applies

import json
from sage.all import Integer


QUADRATIC_POINTS_DATA_PATH = "quadratic_points_catalogue.json"

with open(QUADRATIC_POINTS_DATA_PATH, "r") as qdpts_dat_file:
    qdpts_dat = json.load(qdpts_dat_file)


def try_najman_trbovic_filter(d, N):

    disc_of_quad_field = Integer(d if d % 4 == 1 else 4 * d)
    ram_primes = disc_of_quad_field.prime_divisors()

    data_this_N = qdpts_dat[str(N)]
    unram_primes = data_this_N["unramified_primes"]
    absurd_intersection = set(unram_primes).intersection(set(ram_primes))
    if len(absurd_intersection) > 0:
        return False
    return True


# The following does the verification


def verify_prop_6pt1():

    vals213 = [26, 30, 35, 39, 40, 48, 50]

    vals438 = [26, 28, 30, 35, 39, 48, 50]

    print("Doing verification for Qsqrt213\n")

    for N in vals213:
        N = Integer(N)
        if not try_najman_trbovic_filter(438, N):
            print(f"{N} : najman_trbovic_filter")
        else:
            print(f"{N} is bad")

    print("\nDoing verification for Qsqrt438\n")

    for N in vals438:
        N = Integer(N)
        if not try_najman_trbovic_filter(438, N):
            print(f"{N} : najman_trbovic_filter")
        else:
            print(f"{N} is bad")
