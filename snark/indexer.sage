# -*- coding: utf-8 -*-
"""
Indexer class that takes in an R1CS constraint system and indexes it for Varuna PIOPs
"""

load("algebra/group.sage")
load("algebra/matrix.sage")
load("algebra/field.sage")

class Indexer:

    def __init__(self, A, B, C, z, w=None, x=None):

        z = list(z)
        x = list(x)
        w = list(w)
        A = matrix(A)
        B = matrix(B)
        C = matrix(C)

        # zero pad z and x
        num_constraints_after_padding = self.nearest_power_of_2(matrix(A).nrows())
        num_variables_after_padding = self.nearest_power_of_2(matrix(A).ncols()) # A.ncols() == len(z)
        num_public_inputs_after_padding = self.nearest_power_of_2(len(x))
        z = self.zero_pad_vector_to_length_n(z, num_variables_after_padding)
        #x = self.zero_pad_vector_to_length_n(x, num_variables_after_padding)
        #x = self.zero_pad_vector_to_length_n(x, num_public_inputs_after_padding)

        # pad matrices
        A = zero_pad_matrix_to_n_m(A, num_constraints_after_padding, num_variables_after_padding)
        B = zero_pad_matrix_to_n_m(B, num_constraints_after_padding, num_variables_after_padding)
        C = zero_pad_matrix_to_n_m(C, num_constraints_after_padding, num_variables_after_padding)

        # get variable, constraint, and public input domains
        variable_domain = self.get_subgroup_power_of_2(num_variables_after_padding)
        constraint_domain = self.get_subgroup_power_of_2(num_constraints_after_padding)
        public_input_domain = self.get_subgroup_power_of_2(num_public_inputs_after_padding)
        self.variable_domain = Group(variable_domain)
        self.constraint_domain = Group(constraint_domain)
        self.X = Group(public_input_domain)

        # get non-zero entries subgroups
        n_A = self.nearest_power_of_2(matrix_sparse_norm(A))
        n_B = self.nearest_power_of_2(matrix_sparse_norm(B))
        n_C = self.nearest_power_of_2(matrix_sparse_norm(C))
        K = self.get_subgroup_power_of_2(max(n_A, n_B, n_C))
        K_A = self.get_subgroup_power_of_2(n_A)
        K_B = self.get_subgroup_power_of_2(n_B)
        K_C = self.get_subgroup_power_of_2(n_C)
        self.K = Group(K)
        self.K_A = Group(K_A, ambient=K)
        self.K_B = Group(K_B, ambient=K)
        self.K_C = Group(K_C, ambient=K)

        group_elements_to_file['K'] = self.K.to_list
        group_elements_to_file['K_A'] = self.K_A.to_list
        group_elements_to_file['K_B'] = self.K_B.to_list
        group_elements_to_file['K_C'] = self.K_C.to_list
        group_elements_to_file['C'] = self.variable_domain.to_list
        group_elements_to_file['R'] = self.constraint_domain.to_list
        group_elements_to_file['X'] = self.X.to_list

        self.A = Matrix(A, self.K_A, self.variable_domain, self.constraint_domain, self.X)
        self.B = Matrix(B, self.K_B, self.variable_domain, self.constraint_domain, self.X)
        self.C = Matrix(C, self.K_C, self.variable_domain, self.constraint_domain, self.X)


        matrix_elements_to_file['val_A'] = self.A.val
        matrix_elements_to_file['col_A'] = self.A.col
        matrix_elements_to_file['row_A'] = self.A.row
        matrix_elements_to_file['val_B'] = self.A.val
        matrix_elements_to_file['col_B'] = self.A.col
        matrix_elements_to_file['row_B'] = self.A.row
        matrix_elements_to_file['val_C'] = self.C.val
        matrix_elements_to_file['col_C'] = self.C.col
        matrix_elements_to_file['row_C'] = self.C.row

        # Compute the LDEs for x and w
        x_vec = Vector(x, self.X)
        x_poly = x_vec.low_degree_extension
        x_evals = []
        for h in self.variable_domain.to_list:
            x_evals.append(x_vec.low_degree_extension(x=h))
        w = self.zero_pad_vector_to_length_n(w, len(z) - len(x))
        ratio = len(z) // len(x)
        w_evals = []
        for k in range(0, len(z)):
            if k % ratio == 0:
                w_evals.append(0)
            else:
                w_evals.append(w[k - (k // ratio) - 1] - x_evals[k])
        w_evals = Vector(w_evals, self.variable_domain)
        w_poly = w_evals.low_degree_extension
        w_poly, r = w_poly.quo_rem(self.X.vanishing_polynomial())
        if r!= 0:
            print('Error: Remainder is non-zero.')
            assert(0)

        self.z = Vector(z, self.variable_domain)
        self.x_poly = x_poly
        self.w_poly = w_poly

        #self.x_poly = Vector(x, self.variable_domain).low_degree_extension
        #self.x_poly = Vector(x, self.X).low_degree_extension
        #self.w_poly = self.get_witness_poly(num_variables_after_padding, len(x), w, self.zero_pad_vector_to_length_n(x, num_variables_after_padding))

    def nearest_power_of_2(self, num):
        c = ceil(log(num, 2).n())
        return 2^c

    def zero_pad_vector_to_power_of_2(self, z):
        c = ceil(log(len(z), 2).n())
        nearest_power_of_2 = 2^c
        zeros_on_end = [0] * (nearest_power_of_2 - len(z))
        z = z + zeros_on_end
        return z

    def zero_pad_vector_to_length_n(self, z, n):
        if len(z) > n:
            print('Error: Cannot pad to negative length.')
            assert(0)
        return z + [0] * (n - len(z))

    def get_subgroup_power_of_2(self, n):

        assert(n > 2) # this will return a weird cached object otherwise
        c = log(n, 2).n()
        if c not in ZZ:
            print('Error: n must be a power of 2.')
            assert(0)
        c = ceil(c)
        if n > Fstar.order():
            print('Error: Length of vector is greater than |F*|.')
            assert(0)
        P = get_root_of_unity(c)
        if P == None:
            print('Error: No root of unity of order 2^c exists.')
            assert(0)
        if c > prime_factors[0][1]:
            print('Error: 2^c is not a factor of |F*|.')
            assert(0)
        H = Fstar.subgroup([P])
        return (H)


    def get_witness_poly(self, num_variables, num_public_inputs, w, x):
        w = self.zero_pad_vector_to_length_n(w, num_variables - num_public_inputs)
        ratio = num_variables // num_public_inputs
        w_shifted = []
        for k in range(num_variables):
            if k % ratio == 0:
                w_shifted.append(F(0))
            else:
                w_shifted.append(w[k - (k // ratio) - 1] - x[k])

        w_shifted_poly = Vector(w_shifted, self.variable_domain).low_degree_extension
        w_poly, r = w_shifted_poly.quo_rem(self.X.vanishing_polynomial())

        if r != 0:
            print('Error: Remainder is non-zero.')
            assert(0)

        return w_poly
