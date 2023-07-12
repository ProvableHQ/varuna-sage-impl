from primitives import Group 
from primitives import Vector 
from primitives import Matrix 


class Prover:
    
    #Pre-processing 
    def __init__(self, K, H, A, B, C, x, w):
        
        self.A = Matrix(A)
        self.B = Matrix(B)
        self.C = Matrix(C)
        
        z = x
        z.extend(w)
        self.z = Vector(z) # z = (x, w)
        
        self.K = Group(K)
        self.H = H
        self.K_A = self.A.K
        self.K_B = self.B.K
        self.K_C = self.C.K
        self.z_A = Vector(self.A.to_matrix * self.z.to_vector)
        self.z_B = Vector(self.B.to_matrix * self.z.to_vector)
        self.z_C = Vector(self.C.to_matrix * self.z.to_vector) 