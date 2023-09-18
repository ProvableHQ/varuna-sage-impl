load("algebra/scalar_field.sage")

class Verifier:

    def __init__(self, row_oracles, col_oracles, val_oracles, K, K_A, K_B, K_C, H, X, x, z_poly, x_poly, w_poly):

        self.K = K
        self.K_A = K_A
        self.K_B = K_B
        self.K_C = K_C

        self.H = H
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
        while gamma in self.H.to_list:
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

    def Round_4_rhs(self, sigma, h1, g1, omegas, etas, beta):

        eta_A = etas[0]
        eta_B = etas[1]
        eta_C = etas[2]

        omega_A = omegas[0]
        omega_B = omegas[1]
        omega_C = omegas[2]

        lhs = eta_A * omega_A * self.x_poly(x=beta)
        lhs += eta_A * omega_A * self.X.vanishing_polynomial(x=beta) * self.w_poly(x=beta)
        lhs += eta_B * omega_B * self.x_poly(x=beta)
        lhs += eta_B * omega_B * self.X.vanishing_polynomial(x=beta) * self.w_poly(x=beta)
        lhs += eta_C * omega_C * self.x_poly(x=beta)
        lhs += eta_C * omega_C * self.X.vanishing_polynomial(x=beta) * self.w_poly(x=beta)

        rhs = h1(x=beta) * self.H.vanishing_polynomial(x=beta) + beta * g1(x=beta) + sigma/self.H.order

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
