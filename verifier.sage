
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
            
        eta_A = Fstar.random_element() 
        eta_B = Fstar.random_element() 
        eta_C = Fstar.random_element() 
  
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
        return beta
        
        
