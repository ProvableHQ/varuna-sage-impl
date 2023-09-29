# -*- coding: utf-8 -*-
"""
Field constants for the scalar field used by the BLS-377 elliptic curve
"""

### Prime field modulus for the base field of the BLS-377 elliptic curve group
p = 8444461749428370424248824938781546531375899335154063827935233455917409239041

### Instantiation of GF and its multiplicative subgroup defined by the unit group
F = GF(p)
Fstar = F.unit_group()

### The polynomial ring (with coeffient support of 0..p-1) used to contstruct polynomials used in Varuna proofs
R = PolynomialRing(F, 'x')

### Multiplicative generator of the base field
GEN = F.multiplicative_generator()
Fstar._values = (GEN,)

### The factors of prime field order p-1
prime_factors = factor(p-1)

### 2^47 is a factor of (p-1), thus the TWO_ADICITY defines a multiplicative subgroup of size 2^47
TWO_ADICITY = 47

### The TWO ADIC root of unity which defines a multiplicative generator where (TWO_ADIC_ROOT_OF_UNITY)^(2^TWO_ADICITY) == 1
### This establishes evaluation domains with which to evaluate or interpolate polynomials using Fast-Fourier transforms
ODD_FACTOR = F(60001509534603559531609739528203892656505753216962260608619555)
TWO_ADIC_ROOT_OF_UNITY = Fstar.gen()^ODD_FACTOR

# Returns a primitive '2^n'th root of unity - this helps form the evaluation domain for Varuna polynomials
def get_root_of_unity(n: sage.rings.integer.Integer):
    omega = TWO_ADIC_ROOT_OF_UNITY
    if omega^(2^TWO_ADICITY) != F(1):
        print('Error: Two-adic root of unity is not of order 2^TWO_ADICITY.')
        assert(0)
    return pow(omega, 2^(TWO_ADICITY - n))

# Casts a group G = [g, g^2, ..., g^n=1] into a list of group elements.
def group_to_list(G: sage.groups.group.Group):
    if not isinstance(G, sage.groups.group.Group):
        print('Error: The input G is not a group.')
        assert(0)
    result = []
    g = G.0
    for i in range(0, G.order()):
        result.append(g**i)
    #return result[1:] + [result[0]]
    return result

# Samples a random element from G, which is either of type 'group' or 'subgroup'
# In Sage, there is no method to sample an element from an object of type 'subgroup'
def random_element(G: sage.groups.abelian_gps.abelian_group.AbelianGroup_subgroup):
    if not isinstance(G, sage.groups.group.Group):
        print('Error: The input G is not a group.')
        assert(0)
    g=None
    if isinstance(G, sage.groups.abelian_gps.abelian_group.AbelianGroup_subgroup):
        rand_index = randint(0, G.order() - 1)
        g = group_to_list(G)[rand_index]
    else:
        g = G.random_element()
    return g

# Sorts elements in subgroup of F* by their value.
def sort_by_value(L):
    result = []
    norm_mapping = {}
    for elem in L:
        norm_mapping[F(elem)] = elem
    keys = list(norm_mapping.keys())
    keys.sort()
    for key in keys:
        result.append(norm_mapping[key])
    return result
