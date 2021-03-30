"""utils.py

Some helpful functions
"""

def split_cartan_genus(p):
    """Computes the genus of X_s(p) from Imin Chen's "The Jacobians of
    non-split Cartan modular curves"""

    return (1/24) * (p^2 - 8*p + 11 - 4 * kronecker_symbol(-3,p))

def nonsplit_cartan_genus(p):
    """Computes the genus of X_ns(p) from Imin Chen's "The Jacobians of
    non-split Cartan modular curves"""

    return (1/24) * (p^2 - 10*p + 23 + 6 * kronecker_symbol(-1,p) + 4 * kronecker_symbol(-3,p))


def is_atkin_lehner_divisor(d, N):

    if d.divides(N):
        if gcd(d, N/d) == 1:
            return True
    return False

def c_i_at_2(i, Ndash, N):

    if i == 1:
        if Ndash%4 == 1:
            if is_atkin_lehner_divisor(2, N):
                return 1
            elif 4.divides(N):
                return 0
        if Ndash%4 == 3:
            if is_atkin_lehner_divisor(2,N):
                return 2
            elif 4.divides(N):
                return 3 + kronecker_symbol(-Ndash,2)
            elif 8.divides(N):
                return 3 * (1 + kronecker_symbol(-Ndash,2))

    if i == 2 and Ndash%4 == 3:
        return 1 + kronecker_symbol(-Ndash,2)


def c_i(i, p, Ndash, N):

    if p == 2:
        return c_i_at_2(i, Ndash, N)

    if Ndash%4 == 3:
        return 1 + kronecker_symbol(-Ndash, p)

    else:
        return 1 + kronecker_symbol(-4*Ndash, p)

def fixed_point_number(Ndash,N):

    assert is_atkin_lehner_divisor(Ndash,N), "Ndash does not give an AL involution"

    N_over_Ndash = ZZ(N/Ndash)

    base_contribution = prod([c_i(1, p, Ndash, N) for p in N_over_Ndash.prime_divisors()]) * (-4*Ndash).class_number()

    if Ndash == 2:
        other_contribution = prod([1 + kronecker_symbol(-4,p) for p in ZZ(N/2).prime_divisors()])
    elif Ndash == 3:
        other_contribution = prod([1 + kronecker_symbol(-3,p) for p in ZZ(N/3).prime_divisors()])
    elif Ndash == 4:
        all_divisors_N_over_4 = (ZZ(N/4)).divisors()
        all_AL_divisors = [D for D in all_divisors_N_over_4 if is_atkin_lehner_divisor(D,(ZZ(N/4)))]
        pre_prime_power_divisor_data = [D.is_prime_power(get_data=True) for D in all_AL_divisors]
        prime_power_divisor_data = [(p,nu) for p,nu in pre_prime_power_divisor_data if nu != 0]
        other_contribution = prod([p**(nu // 2) + p**((nu-1) // 2) for p,nu in prime_power_divisor_data])
    elif Ndash%4 == 3:
        other_contribution = prod([c_i(2, p, Ndash, N) for p in N_over_Ndash.prime_divisors()]) * (-Ndash).class_number()
    else:
        other_contribution = 0

    return base_contribution + other_contribution

def genus_of_quotient(N, Ndash):

    return (1/4) * (2 * Gamma0(N).genus() + 2 - fixed_point_number(Ndash, N))


# Some testing of the function, the following examples taken from Galbraith's thesis

assert genus_of_quotient(91,91) == 2
assert genus_of_quotient(46,2) == 3
assert genus_of_quotient(55, 5) == 3
assert genus_of_quotient(56, 2) == 3
assert genus_of_quotient(84,84) == 4
assert genus_of_quotient(92,23) == 1
assert genus_of_quotient(99,99) == 3