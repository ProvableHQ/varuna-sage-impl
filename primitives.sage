p = 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177
F = GF(p)
prime_factors = factor(p-1)
R.<x> = PolynomialRing(F)
Fstar = F.unit_group()


"""
Algebraic Primitives, Part 1 

"""

# Samples a random element from G, which is either of type 'group' or 'subgroup'
# In Sage, there is no method to sample an element from an object of type 'subgroup'
def random_element(G): 
    if not isinstance(G, sage.groups.group.Group): 
        print('Error: The input G is not a group.')
        assert(0)
    if isinstance(G, sage.groups.abelian_gps.abelian_group.AbelianGroup_subgroup): 
        rand_index = randint(0, G.order() - 1)
        g = group_to_list(G)[rand_index]
    else: 
        g = G.random_element()
    return g

# Casts a group G = [g, g^2, ..., g^n=1] into a list of group elements. 
def group_to_list(G):
    if not isinstance(G, sage.groups.group.Group): 
        print('Error: The input G is not a group.')
        assert(0)
    result = []
    g = G.0 
    for i in range(0, G.order()): 
        result.append(g**i)
    return result 

# Returns an element of group G of order 2^r, where r may be composite. 
def element_order_r(G, r): 
    if not isinstance(G, sage.groups.group.Group): 
        print('Error: The input G is not a group.')
        exit(1)
        
    if G.order() % 2^r != 0: 
        print('Error: The order r of the desired element does not divide order(G).')
        exit(1)
        
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
        exit(1)
        
    q,r = f.quo_rem(g) 
    if r!=R(0): 
        print('Error: Remainder should be 0.')
        exit(1)
    if  f != q*g + r: 
        print('Error: Euclidean division failed.')
        exit(1)
    
    return q/q(x=a) 

"""
Algebraic Primitives, Part 2 
Matrix, Vector, and Group Objects.

"""

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
                      

class Group: 
    
    def __init__(self, G): 
        
        if not isinstance(G, sage.groups.group.Group): 
            print('Error: G is not a group object.')
            exit(1)
        
        self.to_group = G
        self.to_list = group_to_list(G)
        self.vanishing_polynomial = vanishing_polynomial(G)
        self.order = G.order()

class Vector: 
                                                                                                  
    def __init__(self, v):  
        
        if not isinstance(v, list): 
            print('Error: v is not a list.')
            exit(1)
        self.to_vector = vector(v)  
        self.H = Group(self.index_group()) 
        
        if len(v) > self.H.order: 
            print('Error: Unable to index as the order of H is less than len(v).')
            exit(1)
        
        self.norm = self.to_vector.norm() #L2 norm
        self.low_degree_extension = self.low_degree_extension()
        
        
    # Returns the index group K of matrix M generated by a kth root of unity, where k is prime. 
    # TODO: In the future, the order of K doesn't have to be prime. We have to find a primitive kth root of unity. 
    # to generate the group. 
    ## NOTE: For testing purposes, the groups in Rust and Sage have to be the same. 
    def index_group(self): 
        
        print('len z', len(self.to_vector))
        n = len(self.to_vector)
        if n > Fstar.order(): 
            print('Error: Length of vector is greater than |F*|.')
            exit(1)
    
        primes = factor(Fstar.order())
        next_prime=None
        
        # set next_prime equal to the smallest prime greater than n. 
        for i in range(len(primes)): 
            if primes[i][0] >= n: 
                print('prime', primes[i][0])
                next_prime = primes[i][0]
                break 
        if next_prime == None: 
            print("Error: No prime found.")
            exit(1)
            
        # find element in F* of order next_prime. 
        H_gen = element_order_p(next_prime, Fstar) 
        H = Fstar.subgroup([H_gen]) # set H = (h)
        return H
    
        
    # Returns the low degree extension of vector v.    
    # For testing purposes, assume there is a prime factor of |F*| which is equal to the length of v 
    def low_degree_extension(self):                                                                                            
        points = []
        for i, h in enumerate(self.H.to_list): 
            points.append((F(h), F(self.to_vector[i])))                                                                                                                                                                         
        f = R.lagrange_polynomial(points)
        return f
