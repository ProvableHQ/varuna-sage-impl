# -*- coding: utf-8 -*-
"""
Prover class which executes the proofs of the Varuna protocol.
"""

load("algebra/field.sage")

class Prover:

    #Pre-processing
    def __init__(self, A, B, C, K, K_A, K_B, K_C, variable_domain, constraint_domain, X, z, x_poly, w_poly):
        self.variable_domain = variable_domain
        self.constraint_domain = constraint_domain
        self.X = X

        self.K = K
        self.K_A = K_A
        self.K_B = K_B
        self.K_C = K_C

        self.A = A
        self.B = B
        self.C = C

        self.x_poly = x_poly
        self.w_poly = w_poly
        self.z_poly = (self.w_poly * self.X.vanishing_polynomial()) + self.x_poly

        self.z = z
        (z_A_prime, z_A, z_B, z_C) = self.z_M()

        self.z_A_lde = z_A
        self.z_B_lde = z_B
        self.z_C_lde = z_C

        output_elements_to_file['z_A_lde'] = self.z_A_lde
        output_elements_to_file['z_B_lde'] = self.z_B_lde
        output_elements_to_file['z_C_lde'] = self.z_C_lde


    def z_M(self):

        z_A_prime = Vector(self.A.to_matrix * self.z.to_vector, self.constraint_domain)
        z_A = Vector(self.A.to_matrix * self.z.to_vector, self.constraint_domain).low_degree_extension
        z_B = Vector(self.B.to_matrix * self.z.to_vector, self.constraint_domain).low_degree_extension
        z_C = Vector(self.C.to_matrix * self.z.to_vector, self.constraint_domain).low_degree_extension

        return (z_A_prime, z_A, z_B, z_C)


    # PIOP 1: Rowcheck
    def Round_1_lhs(self):

        f = self.z_A_lde * self.z_B_lde - self.z_C_lde
        h0, r = f.quo_rem(self.constraint_domain.vanishing_polynomial())
        if r!= 0:
            print('Error: Remainder is non-zero.')
            assert(0)
        if f != h0*self.constraint_domain.vanishing_polynomial() + r:
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
    def Round_3_lhs(self, gamma, etas: list):

        eta_A = etas[0]
        eta_B = etas[1]
        eta_C = etas[2]

        A_z = 0
        B_z = 0
        C_z = 0
        sigma_A = 0
        sigma_B = 0
        sigma_C = 0
        A_z = self.A.bivariate_matrix_polynomial(X=gamma, Y=None)*self.z_poly
        B_z = self.B.bivariate_matrix_polynomial(X=gamma, Y=None)*self.z_poly
        C_z = self.C.bivariate_matrix_polynomial(X=gamma, Y=None)*self.z_poly

        for h in self.variable_domain.to_list:
            sigma_A += A_z(x=h)
            sigma_B += B_z(x=h)
            sigma_C += C_z(x=h)

        sigma = eta_A * sigma_A + eta_B * sigma_B + eta_C * sigma_C
        f = eta_A * A_z + eta_B * B_z + eta_C * C_z
        h_1, r = f.quo_rem(self.variable_domain.vanishing_polynomial()) # h_1 and y * g_1
        if f != h_1*self.variable_domain.vanishing_polynomial() + r:
            print('Error: Division failed')
            assert(0)

        r = r - sigma/self.variable_domain.order

        g_1, s = r.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)])) #divide by y

        if s != 0:
            print('Error: Remainder is non-zero.')
            assert(0)
        if r != g_1*R.lagrange_polynomial([(1, 1), (-1, -1)]) + s:
            print('Error: Division failed')
            assert(0)

        output_elements_to_file['sigma'] = sigma
        output_elements_to_file['h_1'] = h_1
        output_elements_to_file['g_1'] = g_1

        return (sigma, h_1, g_1, sigma_A, sigma_B, sigma_C)

    def Round_4_lhs(self, gamma, beta):

        omega_A = self.A.bivariate_matrix_polynomial(X=gamma)(x=beta) # A_lde(gamma, beta)
        omega_B = self.B.bivariate_matrix_polynomial(X=gamma)(x=beta) # B_lde(gamma, beta)
        omega_C = self.C.bivariate_matrix_polynomial(X=gamma)(x=beta) # C_lde(gamma, beta)

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
        rowcolvalA = self.A.row()*self.A.col()*self.A.val()
        rowcolA = self.A.row()*self.A.col()
        pA = self.constraint_domain.vanishing_polynomial(x=gamma)*self.variable_domain.vanishing_polynomial(x=beta)*self.A.interpolate_on_K(rowcolvalA)
        qA = self.constraint_domain.order*self.variable_domain.order*(gamma*beta - gamma*self.A.col() - beta*self.A.row() + self.A.interpolate_on_K(rowcolA))
        #qA = self.constraint_domain.order*self.variable_domain.order*(gamma - self.A.row())*(beta - self.A.col())
        points_A = []
        for k in self.K_A.to_list:
            points_A.append((F(k), (pA/qA)(x=k)))

        pA_over_qA = R.lagrange_polynomial(points_A)
        hA, sA = (pA - qA*pA_over_qA).quo_rem(self.K_A.vanishing_polynomial)
        xgA = pA_over_qA - omega_A / self.K_A.order
        gA, rA = xgA.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)]))

        if rA != 0:
            print('Error: Remainder rA is not zero.')
            assert(0)
        if pA - qA*pA_over_qA != hA*self.K_A.vanishing_polynomial():
            print('Error')
            assert(0)
        if sA != 0:
            print('Error: Remainder sA is not zero.')
            assert(0)
        if gA.degree() > self.K_A.order or hA.degree() > max(pA.degree(), self.K_A.order - 1 + qA.degree()):
            print('Error: Degree of gA or hA exceeds maximum bound.')
            assert(0)

        ## B
        rowcolvalB = self.B.row()*self.B.col()*self.B.val()
        rowcolB = self.B.row()*self.B.col()
        pB = self.constraint_domain.vanishing_polynomial(x=gamma)*self.variable_domain.vanishing_polynomial(x=beta)*self.B.interpolate_on_K(rowcolvalB)
        qB = self.constraint_domain.order*self.variable_domain.order*(gamma*beta - gamma*self.B.col() - beta*self.B.row() + self.B.interpolate_on_K(rowcolB))
        #qB = self.constraint_domain.order*self.variable_domain.order*(gamma - self.B.row())*(beta - self.B.col())
        points_B = []
        for k in self.K_B.to_list:
            points_B.append((F(k), (pB/qB)(x=k)))

        pB_over_qB = R.lagrange_polynomial(points_B)
        hB, sB = (pB - qB*pB_over_qB).quo_rem(self.K_B.vanishing_polynomial)
        xgB = pB_over_qB - omega_B / self.K_B.order
        gB, rB = xgB.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)]))

        if rB != 0:
            print('Error: Remainder rB is not zero.')
            assert(0)
        if pB - qB*pB_over_qB != hB*self.K_B.vanishing_polynomial():
            print('Error')
            assert(0)
        if sB != 0:
            print('Error: Remainder sB is not zero.')
            assert(0)
        if gB.degree() > self.K_B.order or hB.degree() > max(pB.degree(), self.K_B.order - 1 + qB.degree()):
            print('Error: Degree of gB or hB exceeds maximum bound.')
            assert(0)

        ## C
        rowcolvalC = self.C.row()*self.C.col()*self.C.val()
        rowcolC = self.C.row()*self.C.col()
        pC = self.constraint_domain.vanishing_polynomial(x=gamma)*self.variable_domain.vanishing_polynomial(x=beta)*self.C.interpolate_on_K(rowcolvalC)
        qC = self.constraint_domain.order*self.variable_domain.order*(gamma*beta - gamma*self.C.col() - beta*self.C.row() + self.C.interpolate_on_K(rowcolC))
        #qC = self.constraint_domain.order*self.variable_domain.order*(gamma - self.C.row())*(beta - self.C.col())
        points_C = []
        for k in self.K_C.to_list:
            points_C.append((F(k), (pC/qC)(x=k)))

        pC_over_qC = R.lagrange_polynomial(points_C)
        hC, sC = (pC - qC*pC_over_qC).quo_rem(self.K_C.vanishing_polynomial)
        xgC = pC_over_qC - omega_C / self.K_C.order
        gC, rC = xgC.quo_rem(R.lagrange_polynomial([(1, 1), (-1, -1)]))

        if rC != 0:
            print('Error: Remainder rC is not zero.')
            assert(0)
        if pC - qC*pC_over_qC != hC*self.K_C.vanishing_polynomial():
            print('Error')
            assert(0)
        if sC != 0:
            print('Error: Remainder sC is not zero.')
            assert(0)
        if gC.degree() > self.K_C.order or hC.degree() > max(pC.degree(), self.K_C.order - 1 + qC.degree()):
            print('Error: Degree of gC or hC exceeds maximum bound.')
            assert(0)

        output_elements_to_file['h_a'] = hA
        output_elements_to_file['h_b'] = hB
        output_elements_to_file['h_c'] = hC
        output_elements_to_file['g_a'] = gA
        output_elements_to_file['g_b'] = gB
        output_elements_to_file['g_c'] = gC

        return (hA, hB, hC, gA, gB, gC)

    def Round_6_lhs(self, hs, deltas):
        hA = hs[0]
        hB = hs[1]
        hC = hs[2]

        delta_A = deltas[0]
        delta_B = deltas[1]
        delta_C = deltas[2]

        h2 = delta_A * hA * self.K_A.order/self.K.order
        h2 += delta_B * hB * self.K_B.order/self.K.order
        h2 += delta_C * hC * self.K_C.order/self.K.order

        output_elements_to_file['h_2'] = h2

        return h2
