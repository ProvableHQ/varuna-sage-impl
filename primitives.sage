p = 79
F = GF(p)
R = PolynomialRing(F, 'x')
Fstar = F.unit_group()

"""
Algebraic Primitives, Part 1 

"""

# Returns element an element of group G of order n, where n is prime 
def element_order_p(p: int, G): 
    
    if not isinstance(G, sage.groups.group.Group): 
        print('Error: The input G is not a group.')
        exit(1)
        
    if not p.is_prime(): 
        print('Error: The input p is not prime.')
        exit(1)
    
    if G.order() % p != 0: 
        print('Error: The order p of the desired element does not divide order(G).')
        exit(1)
        
    # finds element of order p 
    g = G.random_element()
    h = g^(G.order()/p) 
    while h.order()==1: 
        g = G.random_element()
        h = g^(G.order()/p) 
    return h 

# Casts a group G = [g, g^2, ..., g^n=1] into a list of ambient group elements. 
def group_to_list(G):
    
    if not isinstance(G, sage.groups.group.Group): 
        print('Error: The input G is not a group.')
        exit(1)
        
    # returns a list out of the group elements 
    result = []
    g = G.0 
    for i in range(0, G.order()): 
        result.append(g**i)
    return result 

#Return the minimal degree vanishing polynomial over the subgroup H of F*. 
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

# Returns the Lagrange polynomial defined over the set S contained in F* at point a \in F*. 
def lagrange_polynomial(S, a): 
    
    if isinstance(S, sage.groups.group.Group): 
        S = group_to_list(S) # 
    if a not in S: 
        print('Error: a is not an element of S.')
        exit(1)
        
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


class Group: 
    
    def __init__(self, G): 
        
        if not isinstance(G, sage.groups.group.Group): 
            print('Error: G is not a group object.')
            exit(1)
        
        self.to_group = G
        self.to_list = group_to_list(G)
        self.vanishing_polynomial = vanishing_polynomial(G)
        self.order = G.order()
             






