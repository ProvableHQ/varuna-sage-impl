# -*- coding: utf-8 -*-
"""
Class to represent A, B, & C in the R1CS constraint system with methods to generate polynomials used in Varuna PIOPs
"""

load("algebra/polynomial.sage")
load("algebra/field.sage")
load("algebra/vector.sage")

class Matrix:

    def __init__(self, M, K, variable_domain, constraint_domain, X):

        self.to_matrix = matrix(M)

        if not isinstance(self.to_matrix, sage.matrix.matrix0.Matrix):
            print('Error: M is not a matrix.')
            assert(0)

        self.K = K
        self.variable_domain = variable_domain
        self.constraint_domain = constraint_domain
        self.X = X
        self.sparse_norm = matrix_sparse_norm(self.to_matrix) # the number of non-zero entries in M

        if matrix_sparse_norm(self.to_matrix) != matrix_sparse_norm(M):
            print('Error: Sparse norm changed after zero padding.')
            assert(0)

        if self.K.order < self.sparse_norm:
            print('Error: The indexing group K is not large enough.')
            assert(0)

        self.K_to_RC = self.K_to_RC()

        self.row = self.row()
        self.col = self.col()
        self.val = self.val()



    # Returns an injective map of elements in K, the indexing group, to a matrix entry index (i, j) \in R x C.
    # Note R = {0, ..., M.nrows}, C = {0, ..., M.ncols}
    def K_to_RC(self):
        mapping = {}
        index=0
        # Maps group elements in K sequentially to a non-zero entry (i,j) in M
        for i in range(self.to_matrix.nrows()):
            for j in range(self.to_matrix.ncols()):
                if self.to_matrix[i,j] > 0:
                    #mapping[self.K.to_list[j]] = (i,j)
                    mapping[self.K.to_list[index]] = (i,j)
                    index+=1

        return mapping


        # Returns the polynomial row_M: K -> H which is constructed by interpolating the points (k, h) \in K x H.
    def row(self):
        points = []
        for k in self.K.to_list:
            h = self.constraint_domain.to_list[0]
            if k in self.K_to_RC.keys():
                h = self.constraint_domain.to_list[self.K_to_RC[k][0]] # maps k to its corresponding row index, then sends the row index to its corresponding element in H
            points.append((F(k), F(h)))

        f = R.lagrange_polynomial(points)
        return f


    # Returns the polynomial col_M: K -> H which is constructed by interpolating the points (k, h) \in K x H.
    def col(self):
        points = []
        for k in self.K.to_list:
            h = self.variable_domain.to_list[0]
            if k in self.K_to_RC.keys():
                col_index = self.K_to_RC[k][1]
                col_index = reindex_by_subdomain(self.variable_domain, self.X, col_index)
                h = self.variable_domain.to_list[col_index] # maps k to its corresponding col index, then sends the col index to its corresponding element in H
            points.append((F(k), F(h)))

        f = R.lagrange_polynomial(points)
        return f



    # Returns the normalized val_M: K -> F polynomial which is constructed by interpolating the points (k, m) and dividing by a constant.
    # In particular, val_M sends k is M[i,j] divided by u_H(row(k), row(k)) * u_H(col(k), col(k))
    def val(self):
        points = []
        for k in self.K.to_list:
            val = 0
            if k in self.K_to_RC.keys(): # do this check because there may be elements of k that don't correspond to any (i,j) since K is zero-padded
                (i, j) = self.K_to_RC[k]
                val = self.to_matrix[i, j]
            points.append((F(k), val))
        f = R.lagrange_polynomial(points)

        return f

    def interpolate_on_K(self, poly):
        points = []
        for k in self.K.to_list:
            points.append((F(k), poly(x=k)))
        return R.lagrange_polynomial(points)


    # Returns the bivariate polynomial representation of matrix M evaluated at either x or y.
    # i.e. this will return a univariate polynomial of the form M(alpha, x) or M(x, beta) or M(alpha, beta).
    def bivariate_matrix_polynomial(self, X=None, Y=None):
        if X == None and Y == None:
            print('Error: X and Y cannot both be None.')
            assert(0)

        num_variables = len(self.variable_domain.to_list)
        num_constraints = len(self.constraint_domain.to_list)
        #transpose = list(matrix(num_variables, num_constraints))

        temp = list(matrix(num_constraints, num_variables))
        #transpose = list(matrix(num_variables, num_constraints))
        for (row_index, row) in enumerate(self.to_matrix):
            for (col_index, val) in enumerate(row):
                c_i = reindex_by_subdomain(self.variable_domain, self.X, col_index)
                #transpose[c_i][row_index] = val
                temp[row_index][c_i] = val


        M_at_alpha_evals = []

        # Return M(alpha, Y)
        if X != None and Y == None:
            for h in self.constraint_domain.to_list:
                M_at_alpha = 0
                for k in self.K_to_RC.keys():
                    (i, j) = self.K_to_RC[k]
                    j = reindex_by_subdomain(self.variable_domain, self.X, j)
                    val = temp[i][j]
                    #val = transpose[i][j]
                    #val = self.to_matrix[i,j]
                    row_at_k = self.row(x=k)
                    col_at_k = self.col(x=k)
                    L_row_at_k = lagrange_basis_polynomial(self.constraint_domain.to_group, row_at_k)(x=X)
                    L_col_at_k = lagrange_basis_polynomial(self.variable_domain.to_group, col_at_k)
                    M_at_alpha += (val*L_row_at_k*L_col_at_k)(x=h)
                M_at_alpha_evals.append((h, M_at_alpha))

        # Return M(X, alpha)
        elif X == None and Y != None:
            for v in self.variable_domain.to_list:
                M_at_alpha = 0
                for k in self.K_to_RC.keys():
                    (i, j) = self.K_to_RC[k]
                    j = reindex_by_subdomain(self.variable_domain, self.X, j)
                    val = temp[i][j]
                    #val = transpose[i][j]
                    #val = self.to_matrix[i,j]
                    row_at_k = self.row(x=k)
                    col_at_k = self.col(x=k)
                    L_row_at_k = lagrange_basis_polynomial(self.constraint_domain.to_group, row_at_k)
                    L_col_at_k = lagrange_basis_polynomial(self.variable_domain.to_group, col_at_k)(x = Y)
                    M_at_alpha += (val*L_row_at_k*L_col_at_k)(x=v)
                M_at_alpha_evals.append((v,M_at_alpha))

        else:
            print("Unsupported option.")
            assert(0)

        f = R.lagrange_polynomial(M_at_alpha_evals)

        return f

#Returns the sparse norm of the matrix M.
def matrix_sparse_norm(M):
    n = 0
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            if M[i, j] > 0:
                n += 1
    return n

def zero_pad_matrix_to_n_m(M, m, n):
    # Pad rows
    for _ in range(m - M.nrows()):
        M = M.insert_row(M.nrows(), [0] * M.ncols())

    # Pad columns
    M = M.transpose()
    for _ in range(n - M.nrows()):
        M = M.insert_row(M.nrows(), [0] * M.ncols())

    return M.transpose()

#Returns the sparse norm of the matrix M.
def matrix_sparse_norm(M):
    n = 0
    for i in range(M.nrows()):
        for j in range(M.ncols()):
            if M[i, j] > 0:
                n += 1
    return n
