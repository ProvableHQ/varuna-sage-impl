# -*- coding: utf-8 -*-
"""
Group type that wraps the SageMath group type to include extra information about the vanishing polynomial and selector
polynomials of the circuit(s) being proved over the group.
"""

load("algebra/polynomial.sage");
load("algebra/field.sage");

class Group:

    def __init__(self, G, gen=None, ambient=None):

        if not isinstance(G, sage.groups.group.Group):
            print('Error: G is not a group object.')
            assert(0)

        self.to_group = G
        self.gen = gen
        self.to_list = group_to_list(G)
        self.vanishing_polynomial = vanishing_polynomial(G)

        # ambient is the ambient group of G, i.e. G <= 'ambient' as groups.
        self.ambient = ambient
        self.ambient_vanishing_polynomial = None
        self.selector = None

        self.order = G.order()

        if ambient != None:
            self.ambient_vanishing_polynomial = vanishing_polynomial(ambient)
            f = self.ambient_vanishing_polynomial
            g = self.vanishing_polynomial
            q,r = f.quo_rem(g)
            if r != R(0):
                print('Error: Remainder is non-zero.')
                assert(0)
            if f != q*g + r:
                print('Error: Division failed.')
                assert(0)
            self.selector = F(G.order() / ambient.order())*q # the selector polynomial

#Reindex constraint matrix columns by the subdomain their oracle polynomials will be defined over.
def reindex_by_subdomain(self: Group, other: Group, index: int):
    period = self.order / other.order
    if index < other.order:
        return index * period
    else:
        i = index - other.order
        x = period - 1
        return i + (i / x) + 1
