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