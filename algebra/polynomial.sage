# -*- coding: utf-8 -*-
"""
Helper functions used to create polynomials used in Varuna PIOPs
"""

#Return the minimal degree vanishing polynomial over the subgroup H of F*.
# v_S = (x-s1)(x-s2)...(x-sn)
def vanishing_polynomial(S):
    if isinstance(S, sage.groups.group.Group):
        S = group_to_list(S)
    prod = R(1)
    for s in S:
        if s == F(0):
            prod *= R.lagrange_polynomial([(1, 1), (-1, -1)]) # = X
        else:
            prod *= R.lagrange_polynomial([(F(s), 0), (0, -F(s))])  # = prod*(X - s)
    return prod

# Return the polynomial v_S / (X - a).
def d_vanishing_polynomial(S, a):
    if isinstance(S, sage.groups.group.Group):
        S = group_to_list(S)
    if a not in S:
        print('Error: a is not an element of S.')
        assert(0)
    v_S = vanishing_polynomial(S)
    x_minus_a = R.lagrange_polynomial([(F(a), 0), (0, -F(a))]) # (X - a)
    q, r = v_S.quo_rem(x_minus_a) # q = v_S / (x-a)
    if r != R(0):
        print('Error: Remainder is not zero.')
        assert(0)
    return q


def lagrange_basis_polynomial(S, a):
    if isinstance(S, sage.groups.group.Group):
        S = group_to_list(S) #
    if a not in S:
        print('Error: a is not an element of S.')
        assert(0)

    f = vanishing_polynomial(S)
    g=None
    if a == F(0):
        g = R.lagrange_polynomial([(1, 1), (-1, -1)]) # = X
    else:
        g = R.lagrange_polynomial([(F(a), 0), (0, -F(a))]) # = (X-a)
    if g == None:
        print('Error: g is None')
        assert(0)

    q,r = f.quo_rem(g)
    if r!=R(0):
        print('Error: Remainder should be 0.')
        assert(0)
    if  f != q*g + r:
        print('Error: Euclidean division failed.')
        assert(0)

    return q*(F(a) / F(len(S)))

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

def interpolate(poly, domain):
    points = []
    for k in domain:
        points.append((F(k), poly(x=k)))
    return R.lagrange_polynomial(points)
