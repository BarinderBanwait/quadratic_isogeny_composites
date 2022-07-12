"""hyperelliptic_verifs.py

This verifies the claims in Section 6 of the paper.

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

# We first show that, for all N as in the statement of Proposition 6.1,
# apart from N = 37, the Trbovic filter method applies

import json
from sage.all import Integer


QUADRATIC_POINTS_DATA_PATH = "quadratic_points_catalogue.json"

with open(QUADRATIC_POINTS_DATA_PATH, "r") as qdpts_dat_file:
    qdpts_dat = json.load(qdpts_dat_file)


def try_trbovic_filter(d, N):

    disc_of_quad_field = Integer(d if d % 4 == 1 else 4 * d)
    ram_primes = disc_of_quad_field.prime_divisors()

    data_this_N = qdpts_dat[str(N)]
    unram_primes = data_this_N["unramified_primes"]
    absurd_intersection = set(unram_primes).intersection(set(ram_primes))
    if len(absurd_intersection) > 0:
        return False
    return True


# The following does the verification


def verify_prop_5pt2():

    vals213 = [26, 30, 35, 39, 40, 48, 50]

    print("Doing verification for Qsqrt213\n")

    for N in vals213:
        N = Integer(N)
        if not try_trbovic_filter(438, N):
            print(f"{N} : trbovic_filter")
        else:
            print(f"{N} is bad")
