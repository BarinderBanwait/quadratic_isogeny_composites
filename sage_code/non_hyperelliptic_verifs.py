"""non_hyperelliptic_verifs.py

    This has the sage code to verify some of the cases in Section 7 of
    the paper.

"""

### N = 43

# We first show that J0(43)_(K) = J0(43)_(Q). This is achieved with the
# following functions from `Isogeny Primes`.


def is_torsion_same(p, chi, B=30, uniform=False):
    """Returns true if the minus part of J0(p) does not gain new torsion when
    base changing to K"""
    M = ModularSymbols(p)
    S = M.cuspidal_subspace()
    T = S.atkin_lehner_operator()
    S_min = (T + parent(T)(1)).kernel()
    J0_min = S_min.abelian_variety()

    d = 2

    if uniform:
        frob_poly_data = [(q, d) for q in prime_range(d + 2, B) if q != p]
    else:
        frob_poly_data = [
            (q, 1) if chi(q) == 1 else (q, d)
            for q in prime_range(d + 2, B)
            if gcd(q, p) == 1
        ]

    point_counts = []

    for q, i in frob_poly_data:
        frob_pol_q = J0_min.frobenius_polynomial(q)
        frob_mat = companion_matrix(frob_pol_q)
        point_counts.append((frob_mat ** i).charpoly()(1))

    # Recall that the rational torsion on J0(p) is entirely contained in
    # the minus part (theorem of Mazur), so checking no-growth of torsion
    # in minus part is done simply as follows

    return J0(p).rational_torsion_order(proof=False) == gcd(point_counts)


def is_rank_of_twist_zero(p, chi):
    """Returns true if the rank of the twist of the minus part of X_0(p)
    by the character chi is zero"""
    ML = ModularSymbols(p, base_ring=chi.base_ring())
    SL = ML.cuspidal_subspace()
    TL = SL.atkin_lehner_operator()
    S_min_L = (TL + parent(TL)(1)).kernel()

    for S in S_min_L.decomposition():
        my_map = S.rational_period_mapping()
        w = ML([0, oo])
        wmap = my_map(w)
        if wmap == 0:
            return False
        tw = ML.twisted_winding_element(0, chi)
        twmap = my_map(tw)
        if twmap == 0:
            return False

    return True


def check_mwgp_same_minus(p, d):
    chi = kronecker_character(d)
    if is_rank_of_twist_zero(p, chi):
        if is_torsion_same(p, chi) == 1:
            return True
    return False


# The following then tests that for p = 43 and d = 213, the MW groups
# are the same

p = 43
d = 213
check_mwgp_same_minus(p, d)

# The following code shows that the only elliptic curve defined over K
# with CM and admitting a 43-isogeny is the j-invariant written in Section
# 7.1

R.<x> = QQ[]
f = x^2 - 213
K.<a> = NumberField(f)
cm_jinvs = cm_j_invariants(K)
output = []
for j in cm_jinvs:
    Ej = EllipticCurve_from_j(j)
    ipd = Ej.isogenies_prime_degree(minimal_models=False)
    pdi = [phi.degree() for phi in ipd]
    if 43 in pdi:
        output.append(j)
output

### N = 65

# Here is the implementation of the Özman sieve from `Isogeny Primes`

def oezman_sieve(p, N):
    """If p is unramified in Q(sqrt(-N)) this always returns True.
    Otherwise returns True iff p is in S_N or . Only makes sense if p ramifies in K"""

    M = QuadraticField(-N)
    if p.divides(M.discriminant()):
        return True

    pp = (M * p).factor()[0][0]
    C_M = M.class_group()
    if C_M(pp).order() == 1:
        return True

    return False

# Here's a wrapper which gets the ps we are able to try

def try_oezman_sieve(d, N):

    if not N.is_squarefree():
        return True

    disc_of_quad_field = d if d%4 == 1 else 4*d
    ram_primes = disc_of_quad_field.prime_divisors()

    for p in ram_primes:
        if not oezman_sieve(p,N):
            return False
    return True

# Now we try it

try_oezman_sieve(213, 65)  # False, which is what we wanted