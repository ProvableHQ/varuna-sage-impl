# -*- coding: utf-8 -*-
"""
Vector class containing the evaluation domain, indexing group, and low degree polyonomial extension representation of
a vector. This is used to represent the witness Vector and the representations of its multiplications with the matrices
A, B, and C representing the constraints of the circuit z_A, z_B, and z_C.
"""

load("algebra/field.sage")
load("algebra/polynomial.sage")

class Vector:

    def __init__(self, v, H):

        self.to_vector = vector(v)
        self.H = H

        if len(v) != self.H.order:
            print('Error: The length of the vector and its indexing group are different.')
            assert(0)

        # self.norm = self.to_vector.norm() #L2 norm
        self.len = len(self.to_vector)
        self.low_degree_extension = self.low_degree_extension()


    # Returns the low degree extension polynomial of vector v.
    # This is done by returning the Lagrange interpolation of the points (h^i, v[i]) for 0 <= i < len(v).
    def low_degree_extension(self):
        points = []
        for i, h in enumerate(self.H.to_list):
            points.append((F(h), F(self.to_vector[i])))
        f = R.lagrange_polynomial(points)
        return f
