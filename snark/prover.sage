load("algebra/scalar_field.sage")

class Prover:

    #Pre-processing
    def __init__(self, A, B, C, K, K_A, K_B, K_C, H, z, W, w_poly, X, x_poly):
        self.H = H
        self.X = X
        self.W = W
        self.K = K
        self.K_A = K_A
        self.K_B = K_B
        self.K_C = K_C

        self.A = A
        self.B = B
        self.C = C
        self.z = z
        self.w_poly = w_poly
        self.x_poly = x_poly
        self.z_poly = (self.w_poly * self.X.vanishing_polynomial()) + self.x_poly

        self.z_A_lde = Vector(self.z_M(A, z), H).low_degree_extension
        self.z_B_lde = Vector(self.z_M(B, z), H).low_degree_extension
        self.z_C_lde = Vector(self.z_M(C, z), H).low_degree_extension

        if self.K_A != self.A.K:
            print('Error: K_A is not the same as A.K.')
            assert(0)
        if self.K_B != self.B.K:
            print('Error: K_B is not the same as B.K.')
            assert(0)
        if self.K_C != self.C.K:
            print('Error: K_C is not the same as C.K.')
            assert(0)

    @staticmethod
    def z_M(M, z):
        z_m = []
        for i in range(0, M.to_matrix.nrows()):
            row = M.to_matrix[i]
            acc = 0
            for j, el in enumerate(row):
                acc += z.to_vector[j]*el
            z_m.append(acc)
        return z_m

    # PIOP 1: Rowcheck
    def Round_1_lhs(self):

        f = self.z_A_lde * self.z_B_lde - self.z_C_lde
        h0, r = f.quo_rem(self.H.vanishing_polynomial())
        if r!= 0:
            print('Error: Remainder is non-zero.')
            assert(0)
        if f != h0*self.H.vanishing_polynomial() + r:
            print('Error: Division failed.')
            assert(0)

        output_elements_to_file['x_lde'] = self.x_poly
        output_elements_to_file['w_lde'] = self.w_poly
        output_elements_to_file['z_lde'] = self.z_poly
        output_elements_to_file['h_0'] = h0

        return h0

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

        M_a = self.A.bivariate_matrix_polynomial(X=gamma, label="a", input_domain=self.X)
        M_b = self.B.bivariate_matrix_polynomial(X=gamma, label="b", input_domain=self.X)
        M_c = self.C.bivariate_matrix_polynomial(X=gamma, label="c", input_domain=self.X)
        z_a = M_a * self.z_poly
        z_b = M_b * self.z_poly
        z_c = M_c * self.z_poly
        sigma_A = 0
        for h in self.H.to_list:
            sigma_A += z_a(h)
        sigma_B = 0
        for h in self.H.to_list:
            sigma_B += z_b(h)
        sigma_C = 0
        for h in self.H.to_list:
            sigma_C += z_c(h)
        sigma = eta_A * sigma_A + eta_B * sigma_B + eta_C * sigma_C

        f = (eta_A * M_a * self.z_poly
             + eta_B * M_b * self.z_poly
             + eta_C * M_c * self.z_poly)

        h_1, r = f.quo_rem(self.H.vanishing_polynomial()) # h_1 and y * g_1

        if f != h_1*self.H.vanishing_polynomial() + r:
            print('Error: Division failed')
            assert(0)

        r = r - sigma/self.H.order

        g_1, s = r.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)]))

        if s != 0:
            print('Error: Remainder is non-zero.')
            assert(0)
        if r != g_1*R.lagrange_polynomial([(1, 1), (-1, -1)]) + s:
            print('Error: Division failed')
            assert(0)

        output_elements_to_file['sigma'] = sigma
        output_elements_to_file['h1'] = h_1
        output_elements_to_file['g1'] = g_1

        return (sigma, h_1, g_1, sigma_A, sigma_B, sigma_C)

    def Round_4_lhs(self, gamma, beta):

        omega_A = self.A.bivariate_matrix_polynomial(X=gamma, label="a", input_domain=self.X)(x=beta)
        omega_B = self.B.bivariate_matrix_polynomial(X=gamma, label="b", input_domain=self.X)(x=beta)
        omega_C = self.C.bivariate_matrix_polynomial(X=gamma, label="c", input_domain=self.X)(x=beta)

        output_elements_to_file['omegaA'] = omega_A
        output_elements_to_file['omegaB'] = omega_B
        output_elements_to_file['omegaC'] = omega_C

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
            print('Error: Remainder rA is not zero.')
            assert(0)
        hA, sA = (pA - qA*fA).quo_rem(self.K_A.vanishing_polynomial())
        if pA - qA*fA != hA*self.K_A.vanishing_polynomial():
            print('Error')
            assert(0)
        if sA != R(0):
            print('Error: Remainder sA is not zero.')
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
            print('Error: Remainder rB is not zero.')
            assert(0)
        hB, sB = (pB - qB*fB).quo_rem(self.K_B.vanishing_polynomial())
        if pB - qB*fB != hB*self.K_B.vanishing_polynomial():
            print('Error: Division failed.')
            assert(0)
        if sB != R(0):
            print('Error: Remainder sB is not zero.')
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
            print('Error: Remainder rC is not zero.')
            assert(0)
        hC, sC = (pC - qC*fC).quo_rem(self.K_C.vanishing_polynomial())
        if pC - qC*fC != hC*self.K_C.vanishing_polynomial():
            print('Error: Division failed.')
            assert(0)
        if sC != R(0):
            print('Error: Remainder sC is not zero.')
            assert(0)
        if gC.degree() > self.K_C.order or hC.degree() > max(pC.degree(), self.K_C.order - 1 + qC.degree()):
            print('Error: Degree of gC or hC exceeds maximum bound.')
            assert(0)

        output_elements_to_file['hA'] = hA
        output_elements_to_file['hB'] = hB
        output_elements_to_file['hC'] = hC
        output_elements_to_file['gA'] = gA
        output_elements_to_file['gB'] = gB
        output_elements_to_file['gC'] = gC

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

        output_elements_to_file['h2'] = h2

        return h2