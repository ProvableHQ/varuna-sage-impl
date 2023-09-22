p = 8444461749428370424248824938781546531375899335154063827935233455917409239041
F = GF(p)
prime_factors = factor(p-1)
R = PolynomialRing(F, 'x')
Fstar = F.unit_group()
GEN = F.multiplicative_generator()
Fstar._values = (GEN,)
TWO_ADICITY = 47
ODD_FACTOR = F(60001509534603559531609739528203892656505753216962260608619555)
TWO_ADIC_ROOT_OF_UNITY = Fstar.gen()^ODD_FACTOR

if ODD_FACTOR * 2^TWO_ADICITY != p-1: 
    print('Error: Two-adicity is incorrect.')
    assert(0)


"""
Algebraic Primitives, Part 1 

"""

randomness_to_file = {}
output_elements_to_file = {}
group_elements_to_file = {} 
matrix_elements_to_file = {}

# Samples a random element from G, which is either of type 'group' or 'subgroup'
# In Sage, there is no method to sample an element from an object of type 'subgroup'
def random_element(G): 
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

"""
def zero_pad_matrix_to_n_m(M, n, m): 
    for _ in range(n - M.nrows()): 
        M = M.insert_row(M.nrows() - 1, [0] * M.ncols())

    M = M.transpose() 
    for _ in range(m - M.nrows()): 
        M = M.insert_row(M.nrows() - 1, [0] * M.ncols())

    return M.transpose()

"""

def zero_pad_matrix_to_n_m(M, m, n):
    # Pad rows
    for _ in range(m - M.nrows()): 
        M = M.insert_row(M.nrows(), [0] * M.ncols())

    # Pad columns
    M = M.transpose()
    for _ in range(n - M.nrows()): 
        M = M.insert_row(M.nrows(), [0] * M.ncols())
    
    return M.transpose()


# Sorts elements in subgroup of F* by their value. 
def sort_by_value(L: list):
    result = []
    norm_mapping = {}
    for elem in L: 
        norm_mapping[F(elem)] = elem
    keys = list(norm_mapping.keys())
    keys.sort() 
    for key in keys: 
        result.append(norm_mapping[key])
    return result 
        
# Casts a group G = [g, g^2, ..., g^n=1] into a list of group elements. 
def group_to_list(G):
    if not isinstance(G, sage.groups.group.Group): 
        print('Error: The input G is not a group.')
        assert(0)
    result = []
    g = G.0
    for i in range(0, G.order()): 
        result.append(g**i)
    #return result[1:] + [result[0]]
    return result 

# Returns an element of group G of order 2^r, where r may be composite. 
def element_order_r(G, r): 
    if not isinstance(G, sage.groups.group.Group): 
        print('Error: The input G is not a group.')
        assert(0)
        
    if G.order() % 2^r != 0: 
        print('Error: The order r of the desired element does not divide order(G).')
        assert(0)
        
    g = random_element(G)
    h = g^(G.order()/2^r) 
    tries = 0
    
    #Ends when h is not the identity and is of order 2^r. 
    while h.order()==1 or not is_order_r(h, r): 
        g = random_element(G)
        h = g^(G.order()/2^r)
        tries +=1 
        if tries > 1000: 
            print('Error: Element order 2^r not found.')
            assert(0)
    return h 

# Returns a primitive '2^n'th root of unity. 
def get_root_of_unity(G, n): 
    omega = TWO_ADIC_ROOT_OF_UNITY 
    if omega^(2^TWO_ADICITY) != F(1): 
        print('Error: Two-adic root of unity is not of order 2^TWO_ADICITY.')
        assert(0)
    return pow(omega, 2^(TWO_ADICITY - n)) 

# Checks if the order of a point P is 2^r. 
# Returns False if 2^r' * P = 1 for 0 < r' < r and True otherwise. 
def is_order_r(P,r): 
    for i in range(r): 
        if F(P^(2^i)) == F(1): 
            return False 
    if F(P^(2^r)) != F(1): 
        return False 
    return True 

#Returns the sparse norm of the matrix M. 
def matrix_sparse_norm(M):
    n = 0
    for i in range(M.nrows()): 
        for j in range(M.ncols()): 
            if M[i, j] > 0: 
                n += 1 
    return n

# Returns a matrix R that is M after zero-padding to an n x n matrix. 
def zero_pad_matrix(M, n): 
    if n < max(M.nrows(), M.ncols()): 
        print('Error: M cannot be zero-padded to an n x n matrix.')
        assert(0)
    R = matrix(n)
    for i in range(M.nrows()): 
        for j in range(M.ncols()): 
            R[i,j] = M[i,j]
    return R 

# Returns a vector w that is v after zero-padding to length n. 
def zero_pad_vector(v, n): 
    assert(v != None)
    if n < len(v):
        print('Error: v cannot be zero-padded to length n.')
        assert(0)
    w = zero_vector(n)
    for i in range(len(v)):
        w[i] = v[i]
    return w 

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


"""
# Returns the Lagrange basis polynomial defined over the set S subset F* at point a in F*. 
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
    
    return q/q(x=a)
"""

def reindex_by_subdomain(self, other, index):
    period = self.order / other.order
    if index < other.order:
        return index * period
    else:
        i = index - other.order
        x = period - 1
        return i + (i / x) + 1

"""

Algebraic Primitives, Part 2 
Matrix, Vector, and Group Objects.

"""

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
    

class Group: 
    
    # ambient is the ambient group of G, i.e. G <= 'ambient' as groups. 
    def __init__(self, G, gen=None, ambient=None): 
        
        if not isinstance(G, sage.groups.group.Group): 
            print('Error: G is not a group object.')
            assert(0)
      
        self.to_group = G
        self.gen = gen
        self.to_list = group_to_list(G)
        self.vanishing_polynomial = vanishing_polynomial(G)
        
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
        group_elements_to_file['variable_domain'] = self.variable_domain.to_list
        group_elements_to_file['constraint_domain'] = self.constraint_domain.to_list
        group_elements_to_file['public_input_domain'] = self.X.to_list

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
        P = get_root_of_unity(Fstar, c)
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


class Prover:
    
    #Pre-processing 
    def __init__(self, A, B, C, K, K_A, K_B, K_C, variable_domain, constraint_domain, X, z, x_poly, w_poly):
        self.variable_domain = variable_domain
        self.constraint_domain = constraint_domain
        self.X = X

        self.K = K
        self.K_A = K_A
        self.K_B = K_B
        self.K_C = K_C

        self.A = A
        self.B = B
        self.C = C

        self.x_poly = x_poly
        self.w_poly = w_poly 
        self.z_poly = (self.w_poly * self.X.vanishing_polynomial()) + self.x_poly

        self.z = z
        (z_A_prime, z_A, z_B, z_C) = self.z_M()

        self.z_A_lde = z_A
        self.z_B_lde = z_B
        self.z_C_lde = z_C

        output_elements_to_file['z_A_lde'] = self.z_A_lde
        output_elements_to_file['z_B_lde'] = self.z_B_lde
        output_elements_to_file['z_C_lde'] = self.z_C_lde


    def z_M(self):  

        z_A_prime = Vector(self.A.to_matrix * self.z.to_vector, self.constraint_domain)
        z_A = Vector(self.A.to_matrix * self.z.to_vector, self.constraint_domain).low_degree_extension
        z_B = Vector(self.B.to_matrix * self.z.to_vector, self.constraint_domain).low_degree_extension
        z_C = Vector(self.C.to_matrix * self.z.to_vector, self.constraint_domain).low_degree_extension

        return (z_A_prime, z_A, z_B, z_C)

            
    # PIOP 1: Rowcheck  
    def Round_1_lhs(self): 

        f = self.z_A_lde * self.z_B_lde - self.z_C_lde
        h0, r = f.quo_rem(self.constraint_domain.vanishing_polynomial())
        if r!= 0: 
            print('Error: Remainder is non-zero.')
            assert(0)
        if f != h0*self.constraint_domain.vanishing_polynomial() + r: 
            print('Error: Division failed.')
            assert(0) 

        output_elements_to_file['x_lde'] = self.x_poly
        output_elements_to_file['w_lde'] = self.w_poly
        output_elements_to_file['z_lde'] = self.z_poly
        output_elements_to_file['h_0'] = h0
        
        return h0 
    
    def Round_2_lhs(self, gamma):
        sigma_A = self.z_A_lde(x=gamma)
        sigma_B = self.z_B_lde(x=gamma)
        sigma_C = self.z_C_lde(x=gamma)

        return (sigma_A, sigma_B, sigma_C)
    
    # PIOP 2: Univariate sumcheck
    def Round_3_lhs(self, gamma, etas: list):
            
        eta_A = etas[0]
        eta_B = etas[1]
        eta_C = etas[2]

        A_z = 0
        B_z = 0 
        C_z = 0
        sigma_A = 0
        sigma_B = 0 
        sigma_C = 0
        A_z = self.A.bivariate_matrix_polynomial(X=gamma, Y=None)*self.z_poly
        B_z = self.B.bivariate_matrix_polynomial(X=gamma, Y=None)*self.z_poly
        C_z = self.C.bivariate_matrix_polynomial(X=gamma, Y=None)*self.z_poly

        for h in self.variable_domain.to_list: 
            sigma_A += A_z(x=h)
            sigma_B += B_z(x=h)
            sigma_C += C_z(x=h)

        sigma = eta_A * sigma_A + eta_B * sigma_B + eta_C * sigma_C
        f = eta_A * A_z + eta_B * B_z + eta_C * C_z
        h_1, r = f.quo_rem(self.variable_domain.vanishing_polynomial()) # h_1 and y * g_1 
        if f != h_1*self.variable_domain.vanishing_polynomial() + r: 
            print('Error: Division failed')
            assert(0)

        r = r - sigma/self.variable_domain.order
        
        g_1, s = r.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)])) #divide by y
        
        if s != 0: 
            print('Error: Remainder is non-zero.')
            assert(0)
        if r != g_1*R.lagrange_polynomial([(1, 1), (-1, -1)]) + s: 
            print('Error: Division failed')
            assert(0) 
        
        output_elements_to_file['sigma'] = sigma 
        output_elements_to_file['h1'] = h_1
        output_elements_to_file['g1'] = g_1
        
        return (sigma, h_1, g_1, sigma_A, sigma_B, sigma_C) 
     
    def Round_4_lhs(self, gamma, beta): 
            
        omega_A = self.A.bivariate_matrix_polynomial(X=gamma)(x=beta) # A_lde(gamma, beta)
        omega_B = self.B.bivariate_matrix_polynomial(X=gamma)(x=beta) # B_lde(gamma, beta)
        omega_C = self.C.bivariate_matrix_polynomial(X=gamma)(x=beta) # C_lde(gamma, beta)
        
        output_elements_to_file['omegaA'] = omega_A 
        output_elements_to_file['omegaB'] = omega_B 
        output_elements_to_file['omegaC'] = omega_C  
        
        return (omega_A, omega_B, omega_C)
    
    
    # PIOP 3: Rational sumcheck
    def Round_5_lhs(self, omegas, gamma, beta): 
            
        omega_A = omegas[0]
        omega_B = omegas[1]
        omega_C = omegas[2]
        
        ## A
        rowcolvalA = self.A.row()*self.A.col()*self.A.val()
        pA = self.constraint_domain.vanishing_polynomial(x=gamma)*self.variable_domain.vanishing_polynomial(x=beta)*rowcolvalA
        qA = self.constraint_domain.order*self.variable_domain.order*(gamma - self.A.row())*(beta - self.A.col())
        points_A = [] 
        for k in self.K_A.to_list:
            points_A.append((F(k), (pA/qA)(x=k)))
        
        pA_over_qA = R.lagrange_polynomial(points_A)
        hA, sA = (pA - qA*pA_over_qA).quo_rem(self.K_A.vanishing_polynomial)
        xgA = pA_over_qA - omega_A / self.K_A.order
        gA, rA = xgA.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)]))
      
        if rA != 0:
            print('Error: Remainder rA is not zero.')
            assert(0)
        if pA - qA*pA_over_qA != hA*self.K_A.vanishing_polynomial(): 
            print('Error')
            assert(0)
        if sA != 0:
            print('Error: Remainder sA is not zero.')
            assert(0)
        if gA.degree() > self.K_A.order or hA.degree() > max(pA.degree(), self.K_A.order - 1 + qA.degree()): 
            print('Error: Degree of gA or hA exceeds maximum bound.')
            assert(0)
        
        ## B
        rowcolvalB = self.B.row()*self.B.col()*self.B.val()
        pB = self.constraint_domain.vanishing_polynomial(x=gamma)*self.variable_domain.vanishing_polynomial(x=beta)*rowcolvalB
        qB = self.constraint_domain.order*self.variable_domain.order*(gamma - self.B.row())*(beta - self.B.col())
        points_B = [] 
        for k in self.K_B.to_list:
            points_B.append((F(k), (pB/qB)(x=k)))
        
        pB_over_qB = R.lagrange_polynomial(points_B)
        hB, sB = (pB - qB*pB_over_qB).quo_rem(self.K_B.vanishing_polynomial)
        xgB = pB_over_qB - omega_B / self.K_B.order
        gB, rB = xgB.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)]))

        if rB != 0:
            print('Error: Remainder rB is not zero.')
            assert(0)
        if pB - qB*pB_over_qB != hB*self.K_B.vanishing_polynomial(): 
            print('Error')
            assert(0)
        if sB != 0:
            print('Error: Remainder sB is not zero.')
            assert(0)
        if gB.degree() > self.K_B.order or hB.degree() > max(pB.degree(), self.K_B.order - 1 + qB.degree()): 
            print('Error: Degree of gB or hB exceeds maximum bound.')
            assert(0)
      
        ## C
        rowcolvalC = self.C.row()*self.C.col()*self.C.val()
        pC = self.constraint_domain.vanishing_polynomial(x=gamma)*self.variable_domain.vanishing_polynomial(x=beta)*rowcolvalC
        qC = self.constraint_domain.order*self.variable_domain.order*(gamma - self.C.row())*(beta - self.C.col())
        points_C = [] 
        for k in self.K_C.to_list:
            points_C.append((F(k), (pC/qC)(x=k)))
        
        pC_over_qC = R.lagrange_polynomial(points_C)
        hC, sC = (pC - qC*pC_over_qC).quo_rem(self.K_C.vanishing_polynomial)
        xgC = pC_over_qC - omega_C / self.K_C.order
        gC, rC = xgC.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)]))

        if rC != 0:
            print('Error: Remainder rC is not zero.')
            assert(0)
        if pC - qC*pC_over_qC != hC*self.K_C.vanishing_polynomial(): 
            print('Error')
            assert(0)
        if sC != 0:
            print('Error: Remainder sC is not zero.')
            assert(0)
        if gC.degree() > self.K_C.order or hC.degree() > max(pC.degree(), self.K_C.order - 1 + qC.degree()): 
            print('Error: Degree of gC or hC exceeds maximum bound.')
            assert(0)
        
        output_elements_to_file['hA'] = hA
        output_elements_to_file['hB'] = hB
        output_elements_to_file['hC'] = hC 
        output_elements_to_file['gA'] = gA
        output_elements_to_file['gB'] = gB
        output_elements_to_file['gC'] = gC 
        
        return (hA, hB, hC, gA, gB, gC)
         
    def Round_6_lhs(self, hs, deltas): 
        hA = hs[0]
        hB = hs[1]
        hC = hs[2]
        
        delta_A = deltas[0]
        delta_B = deltas[1]
        delta_C = deltas[2]
        
        h2 = delta_A * self.K_A.selector * hA * self.K_A.order/self.K.order
        h2 += delta_B * self.K_B.selector * hB * self.K_B.order/self.K.order
        h2 += delta_C * self.K_C.selector * hC * self.K_C.order/self.K.order
        
        #h2, r2 = h2.quo_rem(self.K.vanishing_polynomial()) # divide through by v_K 
        
        output_elements_to_file['h2'] = h2 
        
        return h2
        
class Verifier: 
    
    def __init__(self, row_oracles, col_oracles, val_oracles, K, K_A, K_B, K_C, variable_domain, constraint_domain, X, z_poly, x_poly, w_poly):
        
        self.K = K
        self.K_A = K_A
        self.K_B = K_B
        self.K_C = K_C
        
        self.variable_domain = variable_domain
        self.constraint_domain = constraint_domain
        self.X = X
        
        self.row_A = row_oracles[0]
        self.col_A = col_oracles[0]
        self.val_A = val_oracles[0]
        
        self.row_B = row_oracles[1]
        self.col_B = col_oracles[1]
        self.val_B = val_oracles[1]
        
        self.row_C = row_oracles[2]
        self.col_C = col_oracles[2]
        self.val_C = val_oracles[2]

        self.z_poly = z_poly
        self.x_poly = x_poly
        self.w_poly = w_poly
        
    # PIOP 1: Rowcheck 
    def Round_1_rhs(self):    
        gamma = Fstar.random_element()
        while gamma in self.constraint_domain.to_list: 
            gamma = Fstar.random_element()
            
        eta_A = F(1)
        eta_B = Fstar.random_element()
        eta_C = Fstar.random_element()
        
        randomness_to_file['gamma'] = F(gamma)
        randomness_to_file['eta_A'] = F(eta_A)
        randomness_to_file['eta_B'] = F(eta_B)
        randomness_to_file['eta_C'] = F(eta_C)
  
        return (gamma, eta_A, eta_B, eta_C)
        
    
    def Round_2_rhs(self, sigma_A, sigma_B, sigma_C, h, gamma):
        if sigma_A * sigma_B - sigma_C != h(x=gamma) * self.constraint_domain.vanishing_polynomial(x=gamma): 
            print('Error: Rowcheck verification failed.')
            assert(0) 
            
        return 1 
    
    # PIOP 2: Univariate sumcheck
    def Round_3_rhs(self): 
        beta = Fstar.random_element()
        while beta in self.variable_domain.to_list: 
            beta = Fstar.random_element()
            
        randomness_to_file['beta'] = F(beta)
        return beta
    
    def Round_4_rhs(self, sigma, h1, g1, omegas, etas, beta): 
            
        eta_A = etas[0]
        eta_B = etas[1]
        eta_C = etas[2]

        omega_A = omegas[0]
        omega_B = omegas[1]
        omega_C = omegas[2]

        lhs = eta_A * omega_A * self.z_poly(x=beta)
        lhs += eta_B * omega_B * self.z_poly(x=beta)
        lhs += eta_C * omega_C * self.z_poly(x=beta)

        rhs = h1(x=beta) * self.variable_domain.vanishing_polynomial(x=beta) + beta * g1(x=beta) + sigma/self.variable_domain.order
        
        if lhs != rhs: 
            print('Error: Univariate sumcheck verification failed.')
            assert(0)
                
        return 1 
    
    # PIOP 3: Rational sumcheck 
    def Round_5_rhs(self): 
        
        delta_A = F(1)
        delta_B = F.random_element()
        delta_C = F.random_element()
        
        randomness_to_file['delta_A'] = F(delta_A)
        randomness_to_file['delta_B'] = F(delta_B)
        randomness_to_file['delta_C'] = F(delta_C) 
            
        return (delta_A, delta_B, delta_C)
        
    def Round_6_rhs(self, gs, gamma, beta, deltas, omegas, h2): 
        
        zeta = F.random_element()
        randomness_to_file['zeta'] = F(zeta)
        
        g_A = gs[0]
        g_B = gs[1]
        g_C = gs[2]
            
        delta_A = deltas[0]
        delta_B = deltas[1]
        delta_C = deltas[2]
            
        omega_A = omegas[0]
        omega_B = omegas[1]
        omega_C = omegas[2]
        
        rowcolvalA = self.val_A * self.row_A * self.col_A
        constA = self.constraint_domain.order * self.variable_domain.order
        a_A = self.constraint_domain.vanishing_polynomial(x=gamma) * self.variable_domain.vanishing_polynomial(x=beta) * rowcolvalA
        b_A = constA*(gamma - self.row_A)*(beta - self.col_A)
        lhs = delta_A * self.K_A.selector * (a_A - b_A*(R.lagrange_polynomial([(1, 1), (-1, -1)])*g_A + omega_A / self.K_A.order))
        
        rowcolvalB = self.val_B * self.row_B * self.col_B
        constB = self.constraint_domain.order * self.variable_domain.order
        a_B = self.constraint_domain.vanishing_polynomial(x=gamma) * self.variable_domain.vanishing_polynomial(x=beta) * rowcolvalB
        b_B = constB*(gamma - self.row_B)*(beta - self.col_B)
        lhs += delta_B * self.K_B.selector * (a_B - b_B*(R.lagrange_polynomial([(1, 1), (-1, -1)])*g_B + omega_B / self.K_B.order))
        
        rowcolvalC = self.val_C * self.row_C * self.col_C
        constC = self.constraint_domain.order * self.variable_domain.order
        a_C = self.constraint_domain.vanishing_polynomial(x=gamma) * self.variable_domain.vanishing_polynomial(x=beta) * rowcolvalC
        b_C = constC*(gamma - self.row_C)*(beta - self.col_C)
        lhs += delta_C * self.K_C.selector * (a_C - b_C*(R.lagrange_polynomial([(1, 1), (-1, -1)])*g_C + omega_C / self.K_C.order))
        
        #s, w = lhs.quo_rem(self.K.vanishing_polynomial())

        rhs = h2 * self.K.vanishing_polynomial
        
        if lhs(x=zeta) != rhs(x=zeta): 
            print('Error: Rational sumcheck verification failed.')
            assert(0)

        
        return 1
    
def test_cases(A, B, C, z, w=None, x=None): 

    #Indexer(self, A, B, C, z, w=None, x=None)
    #Prover(self, A, B, C, K, K_A, K_B, K_C, variable_domain, constraint_domain, X, x_poly, w_poly):
    #Verifier(self, row_oracles, col_oracles, val_oracles, K, K_A, K_B, K_C, variable_domain, constraint_domain, X, z_poly, x_poly, w_poly):
    
    I = Indexer(matrix(A), matrix(B), matrix(C), vector(z), vector(w), vector(x))
    print()
    P = Prover(I.A, I.B, I.C, I.K, I.K_A, I.K_B, I.K_C, I.variable_domain, I.constraint_domain, I.X, I.z, I.x_poly, I.w_poly)

    row_oracles = [P.A.row, P.B.row, P.C.row]
    col_oracles = [P.A.col, P.B.col, P.C.col]
    val_oracles = [P.A.val, P.B.val, P.C.val]
    
    V = Verifier(row_oracles, col_oracles, val_oracles, I.K, I.K_A, I.K_B, I.K_C, I.variable_domain, I.constraint_domain, I.X, P.z_poly, P.x_poly, P.w_poly)

    # PIOP 1: Rowcheck  
    h = P.Round_1_lhs()
    (gamma, eta_A, eta_B, eta_C) = V.Round_1_rhs()
    (sigA, sigB, sigC) = P.Round_2_lhs(gamma)
    etas = [eta_A, eta_B, eta_C]
    bit_0 = V.Round_2_rhs(sigA, sigB, sigC, h, gamma)
    print('Result of Rowcheck: ', bit_0)
    
    # PIOP 2: Univariate sumcheck 
    (sigma, h1, g1, sigma_A, sigma_B, sigma_C) = P.Round_3_lhs(gamma, etas)
    sigmas = [sigma_A, sigma_B, sigma_C]
    beta = V.Round_3_rhs()
    (omega_A, omega_B, omega_C) = P.Round_4_lhs(gamma, beta)
    omegas = [omega_A, omega_B, omega_C]
    
    bit_1 = V.Round_4_rhs(sigma, h1, g1, omegas, etas, beta)
    print('Result of Univariate sumcheck: ', bit_1)
    
    #PIOP 3: Ratsumcheck 
    (hA, hB, hC, gA, gB, gC) = P.Round_5_lhs(omegas, gamma, beta)
    hs = [hA, hB, hC]
    gs = [gA, gB, gC]
    (deltaA, deltaB, deltaC) = V.Round_5_rhs() 
    deltas = [deltaA, deltaB, deltaC]
    h2 = P.Round_6_lhs(hs, deltas)
    bit_2 = V.Round_6_rhs(gs, gamma, beta, deltas, omegas, h2)
    print('Result of Rational sumcheck: ', bit_2)


# Generates R1CS instances of m x n matrices where z is of the form [1, b^(d+2), b, b^2, ..., b, b, ...] and d is the multiplicative depth of the circuit
# this follows TestCircuit::gen_rand(..)
def gen_r1cs_instance(m, n, b, d):
    # generate the (public and private) witness vector belonging to our TestCircuit
    # the private part
    dummy_vars = n - 3 - d
    w = [2, 4]
    for i in range(0, dummy_vars):
        w.append(w[0])
    # the public part
    x = [1, 8]
    for d in range(3, d+2):
        x.append(w[1] * x[-1])
    z = x + w

    z = vector(z)
    # Generate constraints as per TestCircuit::generate_constraints(...)
    # We initialize the matrix so we can append to it, and cut off the first row later using submatrix
    A = matrix(zero_vector(m))
    B = matrix(zero_vector(m))
    C = matrix(zero_vector(m))

    num_mul_constraints = d - 1 - 1
    # Insert constraints of the form z[1]*z[2]=z[0]
    # TODO: this is hardcoded, make it dynamic based on TestCircuit
    for i in range(0, n - num_mul_constraints):
        a = zero_vector(m)
        a[4] = 1
        A = A.insert_row(A.nrows(), a)
        
        b = zero_vector(m)
        b[5] = 1
        B = B.insert_row(B.nrows(), b)
        
        c = zero_vector(m)
        c[1] = 1
        C = C.insert_row(C.nrows(), c)

    # Insert constraints of the form z[i-1]*z[2]=z[i]
    z_index = 2
    for i in range(0, num_mul_constraints):
        a = zero_vector(m)
        a[z_index - 1] = 1
        A = A.insert_row(A.nrows(), a)

        b = zero_vector(m)
        b[5] = 1
        B = B.insert_row(B.nrows(), b)
          
        c = zero_vector(m)
        c[z_index] = 1
        C = C.insert_row(C.nrows(), c)

        z_index += 1
    
    # Take submatrices using submatrix(i,j,nr,nc), start at entry (i,j), use nr rows, nc cols
    A = A.submatrix(1, 0, n, m)
    B = B.submatrix(1, 0, n, m)
    C = C.submatrix(1, 0, n, m)
    if (A*z).pairwise_product(B*z) != C*z: 
        print('Error: Invalid R1CS instance.')
        assert(0)

    return (A, B, C, z, w, x)


def main(): 
    args = sys.argv[1:]
    
    if args[0] == 'custom': 
        # custom test case  
        A = matrix([[0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 1, 0, 0], [5, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0]])
        B = matrix([[0, 1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]])
        C = matrix([[0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]])
        x = vector([1, 3, 35])
        w = vector([9, 27, 30, 0])
        z = vector([1, 3, 35, 9, 27, 30, 0])
        test_cases(A, B, C, z, w, x)
    else: 
        m = int(args[0])
        n = int(args[1])
        b = int(args[2])
        d = int(args[3])
        (A, B, C, z, w, x) = gen_r1cs_instance(m, n, b, d)
        A = matrix(A)
        B = matrix(B)
        C = matrix(C)
        test_cases(A, B, C, z, w, x)

    # Write r1cs instance to test file 
    with open('r1cs.txt', 'w') as f_r1cs:
        f_r1cs.write('A')
        f_r1cs.write('\n')
        for i in range(0, A.nrows()): 
            for j in range(0, A.ncols()): 
                f_r1cs.write(str(A[i, j]) + ', ')
            f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('B')
        f_r1cs.write('\n')
        for i in range(0, B.nrows()): 
            for j in range(0, B.ncols()): 
                f_r1cs.write(str(B[i, j]) + ', ')
            f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('C')
        f_r1cs.write('\n')
        for i in range(0, C.nrows()): 
            for j in range(0, C.ncols()): 
                f_r1cs.write(str(C[i, j]) + ', ')
            f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('z')
        f_r1cs.write('\n')
        for i in range(0, len(z)): 
            f_r1cs.write(str(z[i]) + ', ')

    f_r1cs.close()

    # Write randomness and test elements to file  
    with open('randomness.txt', 'w') as f_r:
        for key in randomness_to_file: 
            value = randomness_to_file[key]
            f_r.write(str(key) + ' ' + str(value))
            f_r.write('\n')
            f_r.write('\n')
        
    f_r.close()


    with open('outputs.txt', 'w') as f_t:
        for key in output_elements_to_file: 
            value = output_elements_to_file[key]
            f_t.write(str(key) + ' ') 
            for coeff in R(value): 
                f_t.write(str(coeff) + ',')         
            f_t.write('\n')
            f_t.write('\n')
        
    f_t.close()


    with open('groups.txt', 'w') as f_g:
        for key in group_elements_to_file: 
            group = group_elements_to_file[key]
            f_g.write(str(key) + ' ') 
            for g in group: 
                f_g.write(str(F(g))+ ',')         
            f_g.write('\n')
            f_g.write('\n')

    f_g.close()


if __name__ == "__main__":
    main()