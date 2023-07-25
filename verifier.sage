
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
        
        
