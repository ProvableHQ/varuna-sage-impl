
class Verifier: 
    
    def __init__(self, K, x):
        self.K = Group(K)
        self.H = Group(H)
        
        
        def Round_1_rhs_1(self):    
            gamma = Fstar.random_element()
            while gamma in self.H.to_list(): 
                gamma = Fstar.random_element()
            
            eta_A = Fstar.random_element() 
            eta_B = Fstar.random_element() 
            eta_C = Fstar.random_element() 
  
            return (gamma, eta_A, eta_B, eta_C)
        
        
        def Round_1_rhs_2(self, sigma_A, sigma_B, sigma_C, h, gamma):
            if sigma_A * sigma_B - sigma.C != h(x=gamma) * H.vanishing_polynomial(x=gamma): 
                print('Error: Rowcheck verification failed.')
                exit(1) 