from primitives import Group 
from primitives import Vector 
from primitives import Matrix 


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
        h, r = f.quo_rem(self.H.vanishing_polynomial())
        if r!= 0: 
            print('Error: Remainder is non-zero.')
            assert(0)
        if f != h*self.H.vanishing_polynomial() + r: 
            print('Error: Division failed.')
            assert(0) 
            
        return (self.z.low_degree_extension, h) 
    
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
        
        return (sigma, h_1, g_1) 
     