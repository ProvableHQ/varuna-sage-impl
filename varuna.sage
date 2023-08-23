# https://developer.aleo.org/advanced/the_aleo_curves/edwards_bls12#base-field
p = 8444461749428370424248824938781546531375899335154063827935233455917409239041
# https://developer.aleo.org/advanced/the_aleo_curves/bls12-377#base-field
#258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177
F = GF(p)
prime_factors = factor(p-1)
R = PolynomialRing(F, 'x')
Fstar = F.unit_group()

"""
Algebraic Primitives, Part 1 

"""

randomness_to_file = {}
test_elements_to_file = {}

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
        
    sorted_group = sort_by_value(result) 
    if sorted_group[0] != F(1): 
        print('Error: Identity element is not in group.')
        assert(0)
    return sorted_group[1:] + [sorted_group[0]]

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
    
# Returns the Lagrange polynomial defined over the set S \subset F* at point a \in F*. 
def lagrange_polynomial(S, a): 
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

Algebraic Primitives, Part 2 
Matrix, Vector, and Group Objects.

"""

class Matrix:
    
    def __init__(self, M, K, H):
        
        M = matrix(M)
        self.to_matrix = M
        
        if not isinstance(self.to_matrix, sage.matrix.matrix0.Matrix): 
            print('Error: M is not a matrix.')
            assert(0)
            
        self.K = K 
        self.H = H
        self.sparse_norm = matrix_sparse_norm(self.to_matrix)# the number of non-zero entries in M
        
        if matrix_sparse_norm(self.to_matrix) != matrix_sparse_norm(M): 
            print('Error: Sparse norm changed after zero padding.')
            assert(0)
            
        if self.K.order < self.sparse_norm: 
            print('Error: The indexing group K is not large enough.')
            assert(0)
        
        self.K_to_RC = self.K_to_RC() 
        self.R_to_H = self.R_to_H()
        self.C_to_H = self.C_to_H()
        
        self.row = self.row() 
        self.col = self.col() 
        self.val = self.val()
        
        
    def zero_pad(self, M): 
        n = max(M.nrows(), M.ncols())
        R = matrix(n)
        for i in range(M.nrows()): 
            for j in range(M.ncols()): 
                R[i,j] = M[i,j]
        return R 
    
    # Returns the group K generated by an element which is of order the minimal power of 2 which is at
    # least the sparse norm of M.
    def index_group(self): 
        n = self.sparse_norm
        if n > Fstar.order(): 
            print('Error: Sparse norm of M is greater than |F*|.')
            assert(0)
        c = ceil(log(n, 2).n())
        P = element_order_r(Fstar, c)
        while P == None and c <= prime_factors[0][1]: 
            c+=1 
            P = element_order_r(Fstar, c)
        if c > prime_factors[0][1]: 
            print('Error: 2^c is not a factor of |F*|.')
            assert(0)
        K = Fstar.subgroup([P]) # Set K = (P)
        return K             
    
    # Returns an injective map of elements in K, the indexing group, to a matrix entry index (i, j) \in R x C. 
    # Note R = {0, ..., M.nrows}, C = {0, ..., M.ncols}
    def K_to_RC(self): 
        mapping = {}
        index=0
        # Maps group elements in K sequentially to a non-zero entry (i,j) in M 
        for i in range(self.to_matrix.nrows()): 
            for j in range(self.to_matrix.ncols()): 
                if self.to_matrix[i,j] > 0: 
                    mapping[self.K.to_list[index]] = (i,j)
                    index+=1 
        return mapping 
    
    # Returns an injective map from R to H. 
    def R_to_H(self): 
        mapping = {} 
        # Injectively maps the row index in R to its corresponding group element in H
        for i in range(self.to_matrix.nrows()): 
            mapping[i] = self.H.to_list[i]
        return mapping 
    
    # Returns an injective map from C to H. 
    def C_to_H(self): 
        mapping = {} 
        # Injectively maps the row index in C to its corresponding group element in H
        for i in range(self.to_matrix.ncols()): 
            mapping[i] = self.H.to_list[i]
        return mapping 
    
    # Returns the polynomial row_M: K -> H which is constructed by interpolating the points (k, h) \in K x H. 
    def row(self):   
        points = []
        for k in self.K.to_list: 
            h = random_element(self.H.to_group)
            if k in self.K_to_RC.keys(): 
                h = self.R_to_H[self.K_to_RC[k][0]] # maps k to its corresponding row index, then sends the row index to its corresponding element in H
            points.append((F(k), F(h)))
        f = R.lagrange_polynomial(points)
        return f
    
    # Returns the polynomial col_M: K -> H which is constructed by interpolating the points (k, h) \in K x H. 
    def col(self):   
        points = []
        for k in self.K.to_list:
            h = random_element(self.H.to_group)
            if k in self.K_to_RC.keys(): 
                h = self.C_to_H[self.K_to_RC[k][1]] # maps k to its corresponding col index, then sends the col index to its corresponding element in H
            points.append((F(k), F(h)))
        f = R.lagrange_polynomial(points)
        return f
    
    # Returns the normalized val_M: K -> F polynomial which is constructed by interpolating the points (k, m) and dividing by a constant.
    # In particular, val_M sends k is M[i,j] divided by u_H(row(k), row(k)) * u_H(col(k), col(k))
    def val(self):
        points = []
        for k in self.K.to_list: 
            val = 0 
            if k in self.K_to_RC.keys(): 
                (i, j) = self.K_to_RC[k]
                #val = self.to_matrix[i, j] 
                u_row = self.H.order * self.row(x=k)^(self.H.order - 1)
                u_col = self.H.order * self.col(x=k)^(self.H.order - 1)
                val = self.to_matrix[i, j]  / (F(u_row) * F(u_col))
            points.append((F(k), val))                                                                                                                                                                         
        f = R.lagrange_polynomial(points)
        return f
    
        
    # Returns the bivariate polynomial representation of matrix M evaluated at either x or y.
    # i.e. this will return a univariate polynomial of the form M(alpha, x) or M(x, beta) or M(alpha, beta).
    def bivariate_matrix_polynomial(self, X=None, Y=None): 
            
        if X == None and Y == None: 
            print('Error: X and Y cannot both be None.')
            assert(0)
        acc = 0 
        if X != None and Y == None: 
            for k in self.K.to_list: 
                if k not in self.K_to_RC.keys():  
                    continue 
                f = d_vanishing_polynomial(self.H.to_group, self.row(x=k))(x=X)
                g = d_vanishing_polynomial(self.H.to_group, self.col(x=k))
                acc += self.val(x=k)*f*g         
        elif X == None and Y != None: 
            for k in self.K.to_list: 
                if k not in self.K_to_RC.keys():  
                    continue 
                f = d_vanishing_polynomial(self.H.to_group, self.row(x=k))
                g = d_vanishing_polynomial(self.H.to_group, self.col(x=k))(x=Y)
                acc += self.val(x=k)*f*g         
        else: 
            for k in self.K.to_list: 
                if k not in self.K_to_RC.keys():  
                    continue 
                f = d_vanishing_polynomial(self.H.to_group, self.row(x=k))(x=X)
                g = d_vanishing_polynomial(self.H.to_group, self.col(x=k))(x=Y)
                acc += self.val(x=k)*f*g         
        return acc 
                                                                                                                                      
                                                                                                    
class Vector: 
    
    def __init__(self, v, H):  
      
        self.to_vector = vector(v)  
        self.H = H 
        
        if len(v) > self.H.order: 
            print('Error: Unable to index as the order of H is less than len(v).')
            assert(0)
        
        self.norm = self.to_vector.norm() #L2 norm
        self.len = len(self.to_vector)
        self.low_degree_extension = self.low_degree_extension()


    # Returns the indexing group H generated by an element which is of order 
    # at least the minimal power of 2 which is at least len(v). 
    def index_group(self): 
        n = len(self.to_vector)
        if n > Fstar.order(): 
            print('Error: Length of vector is greater than |F*|.')
            assert(0)
        c = ceil(log(n, 2).n())
        P = element_order_r(Fstar, c)
        while P == None and c <= prime_factors[0][1]: 
            c+=1 
            P = element_order_r(Fstar, c)
        if c > prime_factors[0][1]: 
            print('Error: 2^c is not a factor of |F*|.')
            assert(0)
        H = Fstar.subgroup([P])
        return H 
        
    
    # Returns the low degree extension polynomial of vector v.    
    # This is done by returning the Lagrange interpolation of the points (h^i, v[i]) for 0 <= i < len(v).
    def low_degree_extension(self):                                                                                            
        points = []
        for i, h in enumerate(self.H.to_list): 
            if i < len(self.to_vector):
                points.append((F(h), F(self.to_vector[i])))     
        f = R.lagrange_polynomial(points)
        return f
    

class Group: 
    
    # ambient is the ambient group of G, i.e. G <= 'ambient' as groups. 
    def __init__(self, G, ambient=None): 
        
        if not isinstance(G, sage.groups.group.Group): 
            print('Error: G is not a group object.')
            assert(0)
      
        self.to_group = G
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
    
    def __init__(self, A, B, C, z):
        
        (A, B, C, z) = self.zero_padding(A, B, C, z)
        n_A = matrix_sparse_norm(A)
        n_B = matrix_sparse_norm(B)
        n_C = matrix_sparse_norm(C)
        (K, K_A, K_B, K_C) = self.index_group_matrix(n_A, n_B, n_C)
        
        self.K = Group(K)
        self.K_A = Group(K_A, K)
        self.K_B = Group(K_B, K)
        self.K_C = Group(K_C, K)
        self.H = Group(self.index_group_vector(len(z)))
        print("H generator: ", F(self.H.to_list[0]))
        sorted_list = sort_by_value(self.H.to_list)
        for h in sorted_list:
            print("el in H", F(h))

        print("K_A generator: ", F(self.K_A.to_list[0]))
        sorted_list = sort_by_value(self.K_A.to_list)
        for k_a in sorted_list:
            print("el in K_A", F(k_a))
        
        self.A = Matrix(A, self.K_A, self.H)
        self.B = Matrix(B, self.K_B, self.H)
        self.C = Matrix(C, self.K_C, self.H)
        self.z = Vector(z, self.H)
        
        
    # Note that, at the start, we have M.ncols() = len(z) by assumption (for the matrix-vector product to work). 
    # First, we zero-pad z to the closest power of 2 that is at least len(z). Let z' be z after the zero-padding. 
    # Next, we need to make M a square matrix. If M.nrows() <= len(z'), then we make M a square len(z') x len(z') matrix. 
    # Otherwise (i.e. if M.nrows() > len(z')), we zero-pad z' to be of length M.nrows() and make M a square M.nrows() x M.nrows() matrix. 
    # This ensures that H can index into [0, ..., M.nrows() - 1], [0, ..., M.ncols() - 1], and z'. 
    def zero_padding(self, A, B, C, z):
        H = self.index_group_vector(len(z)) 
        if log(H.order(), 2).n() not in ZZ: 
            print('Error: |H| is not a power of 2.')
            assert(0)
            
        z_prime_len = H.order()
        max_nrows = max(A.nrows(), B.nrows(), C.nrows()) 
        n = max(z_prime_len, max_nrows)
        
        A = zero_pad_matrix(A, n)
        B = zero_pad_matrix(B, n)
        C = zero_pad_matrix(C, n)
        z = zero_pad_vector(z, n)
        
        return (A, B, C, z)
            
    # Returns the indexing group K generated by an element which is of order the minimal power of 2 which is at least |M|. 
    @staticmethod
    def index_group_matrix(n_A, n_B, n_C): 
        n = max(n_A, n_B, n_C)
        if n > Fstar.order():
            print('Error: Maximum sparse norm is greater than |F*|.')
            assert(0)
        
        # Find subgroup K < F* generated by element of order 2^c 
        c = ceil(log(n, 2).n())
        if Fstar.order() % 2^c != 0: 
            print('Error: 2^c does not divide |F*|.')
            assert(0)   
        P = element_order_r(Fstar, c)
        while P == None and c < prime_factors[0][1]: 
            c+=1 
            P = element_order_r(Fstar, c)
        if P == None: 
            print('Error: No element found.')
            assert(0)
            
        K = Fstar.subgroup([P]) 
        
        # Find subgroup K_A < K generated by element of order 2^c_A
        c_A = ceil(log(n_A, 2).n())
        if Fstar.order() % 2^c_A != 0: 
            print('Error: 2^c_A does not divide |F*|.')
            assert(0)
        P_A = element_order_r(K, c_A)
        while P_A == None and c_A < c: 
            c_A +=1 
            P_A = element_order_r(K, c_A)
        if P_A == None: 
            P_A = P
        
        # Find subgroup K_B < K generated by element of order 2^c_B
        c_B = ceil(log(n_B, 2).n())
        if Fstar.order() % 2^c_B != 0: 
            print('Error: 2^c_B does not divide |F*|.')
            assert(0)
        P_B = element_order_r(K, c_B)
        while P_B == None and c_B < c: 
            c_B +=1 
            P_B = element_order_r(K_B, c_B)
        if P_B == None: 
            P_B = P
        
        # Find subgroup K_C < K generated by element of order 2^c_C
        c_C = ceil(log(n_C, 2).n())
        if Fstar.order() % 2^c_C != 0: 
            print('Error: 2^c_C does not divide |F*|.')
            assert(0)
        P_C = element_order_r(K, c_C)
        while P_C == None and c_C < c: 
            c_C +=1 
            P_C = element_order_r(K, c_C)
        if P_C == None: 
            P_C = P
        
        K_A = Fstar.subgroup([P_A]) # order 2^c_A
        K_B = Fstar.subgroup([P_B]) # order 2^c_B
        K_C = Fstar.subgroup([P_C]) # order 2^c_C
        
        return (K, K_A, K_B, K_C)
            
        
    # Returns the indexing group H generated by an element which is of order the minimal power of 2 which is at least len(v).
    # If 2^c > len(v), then v is zero padded to size 2^c. 
    @staticmethod 
    def index_group_vector(n): 
        if n > Fstar.order(): 
            print('Error: Length of vector is greater than |F*|.')
            assert(0)
        c = ceil(log(n, 2).n())
        P = element_order_r(Fstar, c)
        while P == None and c <= prime_factors[0][1]: 
            c+=1 
            P = element_order_r(Fstar, c)
        if c > prime_factors[0][1]: 
            print('Error: 2^c is not a factor of |F*|.')
            assert(0)
        H = Fstar.subgroup([P])
        return H
class Prover:
    
    #Pre-processing 
    def __init__(self, A, B, C, K, K_A, K_B, K_C, H, z):
        self.A = A
        self.B = B
        self.C = C
        self.z = z 

        self.H = H
        self.K = K
        self.K_A = K_A
        self.K_B = K_B
        self.K_C = K_C
        
        if self.K_A != self.A.K: 
            print('Error: K_A is not the same as A.K.')
            assert(0)
        if self.K_B != self.B.K: 
            print('Error: K_B is not the same as B.K.')
            assert(0)
        if self.K_C != self.C.K: 
            print('Error: K_C is not the same as C.K.')
            assert(0)
            
        self.z_A_lde = self.z_M(self.A, self.H, self.z)
        self.z_B_lde = self.z_M(self.B, self.H, self.z)
        self.z_C_lde = self.z_M(self.C, self.H, self.z)
    
    #Return the LDE of the matrix-vector product Mz. 
    #M is an instance of class 'Matrix' and z is an instance of class 'Vector' and H is an instance of class 'Group'.
    @staticmethod
    def z_M(M, H, z):
        acc = 0
        for h in H.to_list:
            acc += M.bivariate_matrix_polynomial(None, h)*z.low_degree_extension(x=h)
        return acc
            
    # PIOP 1: Rowcheck  
    def Round_1_lhs(self): 
        
        f = self.z_A_lde * self.z_B_lde - self.z_C_lde

        for coeff in R(self.z_A_lde): 
            print("z_A_lde: ", coeff)
        for coeff in R(self.z_B_lde): 
            print("z_B_lde: ", coeff)
        for coeff in R(self.z_C_lde): 
            print("z_C_lde: ", coeff)

        h0, r = f.quo_rem(self.H.vanishing_polynomial())
        if r!= 0: 
            print('Error: Remainder is non-zero.')
            assert(0)
        if f != h0*self.H.vanishing_polynomial() + r: 
            print('Error: Division failed.')
            assert(0) 
            
        test_elements_to_file['z_lde'] = self.z.low_degree_extension
        test_elements_to_file['h_0'] = h0
        
        return (self.z.low_degree_extension, h0) 
    
    def Round_2_lhs(self, gamma):
        sigma_A = self.z_A_lde(x=gamma)
        sigma_B = self.z_B_lde(x=gamma)
        sigma_C = self.z_C_lde(x=gamma)
            
        return (sigma_A, sigma_B, sigma_C)
    
    # PIOP 2: Univariate sumcheck
    def Round_3_lhs(self, gamma, etas: list, sigmas: list):
            
        eta_A = etas[0]
        eta_B = etas[1]
        eta_C = etas[2]
            
        sigma_A = sigmas[0]
        sigma_B = sigmas[1]
        sigma_C = sigmas[2]
            
        sigma = eta_A * sigma_A + eta_B * sigma_B + eta_C * sigma_C
                
        f = sigma/self.H.order - (eta_A * self.A.bivariate_matrix_polynomial(gamma) 
                                  + eta_B * self.B.bivariate_matrix_polynomial(gamma) 
                                  + eta_C * self.C.bivariate_matrix_polynomial(gamma)) * self.z.low_degree_extension
        
        h_1, r = f.quo_rem(self.H.vanishing_polynomial()) # h_1 and y * g_1 
        
        if f != h_1*self.H.vanishing_polynomial() + r: 
            print('Error: Division failed')
            assert(0) 
            
        g_1, s = r.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)]))
        
        if s != 0: 
            print('Error: Remainder is non-zero.')
            assert(0)
        if r != g_1*R.lagrange_polynomial([(1, 1), (-1, -1)]) + s: 
            print('Error: Division failed')
            assert(0) 
        
        test_elements_to_file['sigma'] = sigma 
        test_elements_to_file['h1'] = h_1
        test_elements_to_file['g1'] = g_1
        
        return (sigma, h_1, g_1) 
     
    def Round_4_lhs(self, gamma, beta): 
            
        omega_A = self.A.bivariate_matrix_polynomial(gamma)(x=beta)
        omega_B = self.B.bivariate_matrix_polynomial(gamma)(x=beta)
        omega_C = self.C.bivariate_matrix_polynomial(gamma)(x=beta)
        
        test_elements_to_file['omegaA'] = omega_A 
        test_elements_to_file['omegaB'] = omega_B 
        test_elements_to_file['omegaC'] = omega_C  
        
        return (omega_A, omega_B, omega_C)
    
    # PIOP 3: Rational sumcheck
    def Round_5_lhs(self, omegas, gamma, beta): 
            
        omega_A = omegas[0]
        omega_B = omegas[1]
        omega_C = omegas[2]
        
        ## A
        pA = self.H.vanishing_polynomial(x=gamma)*self.H.vanishing_polynomial(x=beta)*self.A.val()
        qA = (gamma - self.A.row())*(beta - self.A.col())
        points_A = [] 
        for k in self.K_A.to_list:
            points_A.append((F(k), (pA/qA)(x=k)))
        xgA = R.lagrange_polynomial(points_A) - omega_A / self.K_A.order
        gA, rA = xgA.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)]))
        fA = xgA + omega_A / self.K_A.order
        if rA != R(0):
            print('Error: Remainder is not zero.')
            assert(0)
        hA, sA = (pA - qA*fA).quo_rem(self.K_A.vanishing_polynomial())
        if pA - qA*fA != hA*self.K_A.vanishing_polynomial(): 
            print('Error')
            assert(0)
        if sA != R(0):
            print('Error: Remainder is not zero.')
            assert(0)
        if gA.degree() > self.K_A.order or hA.degree() > max(pA.degree(), self.K_A.order - 1 + qA.degree()): 
            print('Error: Degree of gA or hA exceeds maximum bound.')
            assert(0)
        
        ## B
        pB = self.H.vanishing_polynomial(x=gamma)*self.H.vanishing_polynomial(x=beta)*self.B.val()
        qB = (gamma - self.B.row())*(beta - self.B.col())
        points_B = [] 
        for k in self.K_B.to_list:
            points_B.append((F(k), (pB/qB)(x=k)))
        xgB = R.lagrange_polynomial(points_B) - omega_B / self.K_B.order
        gB, rB = xgB.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)]))
        fB = xgB + omega_B / self.K_B.order
        if rB != R(0):
            print('Error: Remainder is not zero.')
            assert(0) 
        hB, sB = (pB - qB*fB).quo_rem(self.K_B.vanishing_polynomial())
        if pB - qB*fB != hB*self.K_B.vanishing_polynomial(): 
            print('Error: Division failed.')
            assert(0)
        if sB != R(0):
            print('Error: Remainder is not zero.')
            assert(0)
        if gB.degree() > self.K_B.order or hB.degree() > max(pB.degree(), self.K_B.order - 1 + qB.degree()): 
            print('Error: Degree of gB or hB exceeds maximum bound.')
            assert(0)
      
        ## C
        pC = self.H.vanishing_polynomial(x=gamma)*self.H.vanishing_polynomial(x=beta)*self.C.val()
        qC = (gamma - self.C.row())*(beta - self.C.col())
        points_C = [] 
        for k in self.K_C.to_list:
            points_C.append((F(k), (pC/qC)(x=k)))
        xgC = R.lagrange_polynomial(points_C) - omega_C / self.K_C.order
        gC, rC = xgC.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)]))
        fC = xgC + omega_C / self.K_C.order
        if rC != R(0):
            print('Error: Remainder is not zero.')
            assert(0)
        hC, sC = (pC - qC*fC).quo_rem(self.K_C.vanishing_polynomial())
        if pC - qC*fC != hC*self.K_C.vanishing_polynomial(): 
            print('Error: Division failed.')
            assert(0)
        if sC != R(0):
            print('Error: Remainder is not zero.')
            assert(0)
        if gC.degree() > self.K_C.order or hC.degree() > max(pC.degree(), self.K_C.order - 1 + qC.degree()): 
            print('Error: Degree of gC or hC exceeds maximum bound.')
            assert(0)
        
        test_elements_to_file['hA'] = hA
        test_elements_to_file['hB'] = hB
        test_elements_to_file['hC'] = hC 
        test_elements_to_file['gA'] = gA
        test_elements_to_file['gB'] = gB
        test_elements_to_file['gC'] = gC 
        
        return (hA, hB, hC, gA, gB, gC)
        
    # NOTE: I changed this from our spec.  
    def Round_6_lhs(self, hs, deltas): 
        hA = hs[0]
        hB = hs[1]
        hC = hs[2]
        
        delta_A = deltas[0]
        delta_B = deltas[1]
        delta_C = deltas[2]
        
        h2 = delta_A * self.K_A.selector * hA * self.K_A.vanishing_polynomial()
        h2 += delta_B * self.K_B.selector * hB * self.K_B.vanishing_polynomial()
        h2 += delta_C * self.K_C.selector * hC * self.K_C.vanishing_polynomial()
        
        h2, r2 = h2.quo_rem(self.K.vanishing_polynomial()) # divide through by v_K 
        
        test_elements_to_file['h2'] = h2 
        
        return h2
class Verifier: 
    
    def __init__(self, row_oracles, col_oracles, val_oracles, K, K_A, K_B, K_C, H, x):
        
        self.K = K
        self.K_A = K_A
        self.K_B = K_B
        self.K_C = K_C
        
        self.H = H
        
        self.row_A = row_oracles[0]
        self.col_A = col_oracles[0]
        self.val_A = val_oracles[0]
        
        self.row_B = row_oracles[1]
        self.col_B = col_oracles[1]
        self.val_B = val_oracles[1]
        
        self.row_C = row_oracles[2]
        self.col_C = col_oracles[2]
        self.val_C = val_oracles[2]
        
    # PIOP 1: Rowcheck 
    def Round_1_rhs(self):    
        gamma = Fstar.random_element()
        while gamma in self.H.to_list: 
            gamma = Fstar.random_element()
            
        #eta_A = Fstar.random_element() 
        eta_A = F(1)
        eta_B = Fstar.random_element() 
        eta_C = Fstar.random_element() 
        
        randomness_to_file['gamma'] = F(gamma)
        randomness_to_file['eta_A'] = F(eta_A)
        randomness_to_file['eta_B'] = F(eta_B)
        randomness_to_file['eta_C'] = F(eta_C)
  
        return (gamma, eta_A, eta_B, eta_C)
        
    
    def Round_2_rhs(self, sigma_A, sigma_B, sigma_C, h, gamma):
        if sigma_A * sigma_B - sigma_C != h(x=gamma) * self.H.vanishing_polynomial(x=gamma): 
            print('Error: Rowcheck verification failed.')
            assert(0) 
            
        return 1 
    
    # PIOP 2: Univariate sumcheck
    def Round_3_rhs(self): 
        beta = Fstar.random_element()
        while beta in self.H.to_list: 
            beta = Fstar.random_element()
            
        randomness_to_file['beta'] = F(beta)
        return beta
    
    def Round_4_rhs(self, z, sigma, h1, g1, omegas, etas, beta): 
            
        eta_A = etas[0]
        eta_B = etas[1]
        eta_C = etas[2]
            
        omega_A = omegas[0]
        omega_B = omegas[1]
        omega_C = omegas[2]
        
        lhs = sigma/self.H.order - eta_A * omega_A * z(x=beta) 
        lhs -= eta_B * omega_B * z(x=beta)
        lhs -= eta_C * omega_C * z(x=beta)
        
            
        rhs = h1(x=beta) * self.H.vanishing_polynomial(x=beta) + beta * g1(x=beta)
        
        if lhs != rhs: 
            print('Error: Univariate sumcheck verification failed.')
            assert(0)
                
        return 1 
    
    # PIOP 3: Rational sumcheck 
    def Round_5_rhs(self): 
            
        #delta_A = F.random_element()
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
                
        a_A = self.H.vanishing_polynomial(x=gamma) * self.H.vanishing_polynomial(x=beta) * self.val_A
        b_A = (gamma - self.row_A)*(beta - self.col_A)
        lhs = delta_A * self.K_A.selector * (a_A - b_A*(R.lagrange_polynomial([(1, 1), (-1, -1)])*g_A + omega_A / self.K_A.order))
        
        a_B = self.H.vanishing_polynomial(x=gamma) * self.H.vanishing_polynomial(x=beta) * self.val_B
        b_B = (gamma - self.row_B)*(beta - self.col_B)
        lhs += delta_B * self.K_B.selector * (a_B - b_B*(R.lagrange_polynomial([(1, 1), (-1, -1)])*g_B + omega_B / self.K_B.order))
        
        a_C = self.H.vanishing_polynomial(x=gamma) * self.H.vanishing_polynomial(x=beta) * self.val_C
        b_C = (gamma - self.row_C)*(beta - self.col_C)
        lhs += delta_C * self.K_C.selector * (a_C - b_C*(R.lagrange_polynomial([(1, 1), (-1, -1)])*g_C + omega_C / self.K_C.order))
        
        s, w = lhs.quo_rem(self.K.vanishing_polynomial())
        
        rhs = h2 * self.K.vanishing_polynomial
        
        if lhs(x=zeta) != rhs(x=zeta): 
            print('Error: Rational sumcheck verification failed.')
            assert(0)
        
        return 1

def test_cases(A, B, C, z, w=None, x=None): 
    
    I = Indexer(matrix(A), matrix(B), matrix(C), vector(z))
    P = Prover(I.A, I.B, I.C, I.K, I.K_A, I.K_B, I.K_C, I.H, I.z)
    
    row_oracles = [P.A.row, P.B.row, P.C.row]
    col_oracles = [P.A.col, P.B.col, P.C.col]
    val_oracles = [P.A.val, P.B.val, P.C.val]
    
    V = Verifier(row_oracles, col_oracles, val_oracles, I.K, I.K_A, I.K_B, I.K_C, I.H, x)

    # PIOP 1: Rowcheck  
    (zlde, h) = P.Round_1_lhs()
    (gamma, eta_A, eta_B, eta_C) = V.Round_1_rhs()
    (sigA, sigB, sigC) = P.Round_2_lhs(gamma)
    etas = [eta_A, eta_B, eta_C]
    sigmas = [sigA, sigB, sigC]
    bit_0 = V.Round_2_rhs(sigA, sigB, sigC, h, gamma)
    print('Result of Rowcheck: ', bit_0)
    
    # PIOP 2: Univariate sumcheck 
    (sigma, h1, g1) = P.Round_3_lhs(gamma, etas, sigmas)
    beta = V.Round_3_rhs()
    (omega_A, omega_B, omega_C) = P.Round_4_lhs(gamma, beta)
    omegas = [omega_A, omega_B, omega_C]
    
    bit_1 = V.Round_4_rhs(zlde, sigma, h1, g1, omegas, etas, beta)
    print('Result of Univariate sumcheck: ', bit_1)
    
    # PIOP 3: Ratsumcheck 
    (hA, hB, hC, gA, gB, gC) = P.Round_5_lhs(omegas, gamma, beta)
    hs = [hA, hB, hC]
    gs = [gA, gB, gC]
    (deltaA, deltaB, deltaC) = V.Round_5_rhs() 
    deltas = [deltaA, deltaB, deltaC]
    h2 = P.Round_6_lhs(hs, deltas)
    bit_2 = V.Round_6_rhs(gs, gamma, beta, deltas, omegas, h2)
    print('Result of Rational sumcheck: ', bit_2)



# Generates R1CS instances of n x m matrices where z is of the form [1, b^(d+2), b, b^2, ..., b, b, ...] and d is the multiplicative depth of the circuit
def gen_r1cs_instance(n, m, b, d):

    # padded_public_variables: [1, 8]
    # private_variables: [2, 4, 2, 2]

    num_mul_constraints = d - 1

    # first, we generate the (public and private) witness vector belonging to our TestCircuit
    # the public part, consisting of hardcoded 1 and the largest value
    z = [1, b^(d+2)]
    # the private part, consisting of increasing values
    for i in range(1, d+2):    
        z.append(b^i)
    # the private part, consisting of the base value
    for i in range(0, n - num_mul_constraints - 4):    
        z.append(b)

    print("z", z)
    x = z[0:2]
    w = z[2:]
    print("x: ",x)
    print("w: ",w)
    X = Group(Indexer.index_group_vector(len(z)))
    W = Group(Indexer.index_group_vector(len(z)))
    for coeff in R(Vector(vector(x), X).low_degree_extension()): # TODO: use correct field size 
        print("x_lde: ", coeff)
    for coeff in R(Vector(vector(w), W).low_degree_extension()): # TODO: use correct field size 
        print("w_lde: ", coeff)
        
    z = vector(z)#[b^i for i in range(1, m+1)])
    # TODO: clean up snarkVM ordering
    # tmp = z[0]
    # z[0] = z[2]
    # z[2] = tmp
    # We initialize the matrix so we can append to it, and cut off the first row later using submatrix
    A = matrix(zero_vector(m))
    B = matrix(zero_vector(m))
    C = matrix(zero_vector(m))

    num_mul_constraints = d - 1
    # Insert constraints of the form z[1]*z[2]=z[0]
    for i in range(0, n - num_mul_constraints):
        a = zero_vector(m)
        a[2] = 1
        A = A.insert_row(A.nrows(), a)
        
        b = zero_vector(m)
        b[3] = 1
        B = B.insert_row(B.nrows(), b)
        
        c = zero_vector(m)
        c[1] = 1
        C = C.insert_row(C.nrows(), c)

    # Insert constraints of the form z[i-1]*z[2]=z[i]
    z_index = 2
    for i in range(0, num_mul_constraints):
        # TODO: I might need to have a special case for the first mul_depth. Would be nice to clean that weird order up a bit in snarkVM

        a = zero_vector(m)
        if i == 0:
            a[0] = 1
        else:
            a[z_index] = 1
        A = A.insert_row(A.nrows(), a)
        
        b = zero_vector(m)
        b[2] = 1
        B = B.insert_row(B.nrows(), b)
        
        c = zero_vector(m)
        if i == 0:
            c[3] = 1
        else:
            c[z_index + 1] = 1
        C = C.insert_row(C.nrows(), c)
        z_index += 1
    
    # Take submatrices using submatrix(i,j,nr,nc), start at entry (i,j), use nr rows, nc cols
    A = A.submatrix(1, 0, n, m)
    B = B.submatrix(1, 0, n, m)
    C = C.submatrix(1, 0, n, m)

    print("A: ", A)
    print("B: ", B)
    print("C: ", C)

    if (A*z).pairwise_product(B*z) != C*z: 
        print('Error: Invalid R1CS instance.')
        assert(0)

    return (A, B, C, z)


def main(): 
    args = sys.argv[1:]
    n = int(args[0])
    m = int(args[1])
    b = int(args[2])
    d = int(args[3])
    
    (A, B, C, z) = gen_r1cs_instance(n, m, b, d)

    test_cases(A, B, C, z)

    A = matrix(A)
    B = matrix(B)
    C = matrix(C)
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


    with open('test.txt', 'w') as f_t:
        for key in test_elements_to_file: 
            value = test_elements_to_file[key]
            f_t.write(str(key) + ' ') 
            for coeff in R(value): 
                f_t.write(str(coeff) + ',')         
            f_t.write('\n')
            f_t.write('\n')
        
    f_t.close()


if __name__ == "__main__":
    main()




