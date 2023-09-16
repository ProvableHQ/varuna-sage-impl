from sage.all import *
attach("algebra/scalar_field.sage")

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
            assert(_sage_const_0 )
        if self.K_B != self.B.K:
            print('Error: K_B is not the same as B.K.')
            assert(_sage_const_0 )
        if self.K_C != self.C.K:
            print('Error: K_C is not the same as C.K.')
            assert(_sage_const_0 )

        self.z_A_lde = self.z_M(self.A, self.H, self.z)
        self.z_B_lde = self.z_M(self.B, self.H, self.z)
        self.z_C_lde = self.z_M(self.C, self.H, self.z)

    #Return the LDE of the matrix-vector product Mz.
    #M is an instance of class 'Matrix' and z is an instance of class 'Vector' and H is an instance of class 'Group'.
    @staticmethod
    def z_M(M, H, z):
        acc = _sage_const_0
        for h in H.to_list:
            acc += M.bivariate_matrix_polynomial(None, h)*z.low_degree_extension(x=h)
        return acc

        # PIOP 1: Rowcheck
    def Round_1_lhs(self):

        f = self.z_A_lde * self.z_B_lde - self.z_C_lde
        h0, r = f.quo_rem(self.H.vanishing_polynomial())
        if r!= _sage_const_0 :
            print('Error: Remainder is non-zero.')
            assert(_sage_const_0 )
        if f != h0*self.H.vanishing_polynomial() + r:
            print('Error: Division failed.')
            assert(_sage_const_0 )

        output_elements_to_file['z_lde'] = self.z.low_degree_extension
        output_elements_to_file['h_0'] = h0

        return (self.z.low_degree_extension, h0)

    def Round_2_lhs(self, gamma):
        sigma_A = self.z_A_lde(x=gamma)
        sigma_B = self.z_B_lde(x=gamma)
        sigma_C = self.z_C_lde(x=gamma)

        return (sigma_A, sigma_B, sigma_C)

    # PIOP 2: Univariate sumcheck
    def Round_3_lhs(self, gamma, etas: list, sigmas: list):

        eta_A = etas[_sage_const_0 ]
        eta_B = etas[_sage_const_1 ]
        eta_C = etas[_sage_const_2 ]

        sigma_A = sigmas[_sage_const_0 ]
        sigma_B = sigmas[_sage_const_1 ]
        sigma_C = sigmas[_sage_const_2 ]

        sigma = eta_A * sigma_A + eta_B * sigma_B + eta_C * sigma_C

        f = sigma/self.H.order - (eta_A * self.A.bivariate_matrix_polynomial(gamma)
                                  + eta_B * self.B.bivariate_matrix_polynomial(gamma)
                                  + eta_C * self.C.bivariate_matrix_polynomial(gamma)) * self.z.low_degree_extension

        h_1, r = f.quo_rem(self.H.vanishing_polynomial()) # h_1 and y * g_1

        if f != h_1*self.H.vanishing_polynomial() + r:
            print('Error: Division failed')
            assert(_sage_const_0 )

        g_1, s = r.quo_rem(R.lagrange_polynomial([(_sage_const_1 , _sage_const_1 ), (-_sage_const_1 , -_sage_const_1 )]))

        if s != _sage_const_0 :
            print('Error: Remainder is non-zero.')
            assert(_sage_const_0 )
        if r != g_1*R.lagrange_polynomial([(_sage_const_1 , _sage_const_1 ), (-_sage_const_1 , -_sage_const_1 )]) + s:
            print('Error: Division failed')
            assert(_sage_const_0 )

        output_elements_to_file['sigma'] = sigma
        output_elements_to_file['h1'] = h_1
        output_elements_to_file['g1'] = g_1

        return (sigma, h_1, g_1)

    def Round_4_lhs(self, gamma, beta):

        omega_A = self.A.bivariate_matrix_polynomial(gamma)(x=beta)
        omega_B = self.B.bivariate_matrix_polynomial(gamma)(x=beta)
        omega_C = self.C.bivariate_matrix_polynomial(gamma)(x=beta)

        output_elements_to_file['omegaA'] = omega_A
        output_elements_to_file['omegaB'] = omega_B
        output_elements_to_file['omegaC'] = omega_C

        return (omega_A, omega_B, omega_C)

    # PIOP 3: Rational sumcheck
    def Round_5_lhs(self, omegas, gamma, beta):

        omega_A = omegas[_sage_const_0 ]
        omega_B = omegas[_sage_const_1 ]
        omega_C = omegas[_sage_const_2 ]

        ## A
        pA = self.H.vanishing_polynomial(x=gamma)*self.H.vanishing_polynomial(x=beta)*self.A.val()
        qA = (gamma - self.A.row())*(beta - self.A.col())
        points_A = []
        for k in self.K_A.to_list:
            points_A.append((F(k), (pA/qA)(x=k)))
        xgA = R.lagrange_polynomial(points_A) - omega_A / self.K_A.order
        gA, rA = xgA.quo_rem(R.lagrange_polynomial([(_sage_const_1 , _sage_const_1 ), (-_sage_const_1 , -_sage_const_1 )]))
        fA = xgA + omega_A / self.K_A.order
        if rA != R(_sage_const_0 ):
            print('Error: Remainder is not zero.')
            assert(_sage_const_0 )
        hA, sA = (pA - qA*fA).quo_rem(self.K_A.vanishing_polynomial())
        if pA - qA*fA != hA*self.K_A.vanishing_polynomial():
            print('Error')
            assert(_sage_const_0 )
        if sA != R(_sage_const_0 ):
            print('Error: Remainder is not zero.')
            assert(_sage_const_0 )
        if gA.degree() > self.K_A.order or hA.degree() > max(pA.degree(), self.K_A.order - _sage_const_1  + qA.degree()):
            print('Error: Degree of gA or hA exceeds maximum bound.')
            assert(_sage_const_0 )

        ## B
        pB = self.H.vanishing_polynomial(x=gamma)*self.H.vanishing_polynomial(x=beta)*self.B.val()
        qB = (gamma - self.B.row())*(beta - self.B.col())
        points_B = []
        for k in self.K_B.to_list:
            points_B.append((F(k), (pB/qB)(x=k)))
        xgB = R.lagrange_polynomial(points_B) - omega_B / self.K_B.order
        gB, rB = xgB.quo_rem(R.lagrange_polynomial([(_sage_const_1 , _sage_const_1 ), (-_sage_const_1 , -_sage_const_1 )]))
        fB = xgB + omega_B / self.K_B.order
        if rB != R(_sage_const_0 ):
            print('Error: Remainder is not zero.')
            assert(_sage_const_0 )
        hB, sB = (pB - qB*fB).quo_rem(self.K_B.vanishing_polynomial())
        if pB - qB*fB != hB*self.K_B.vanishing_polynomial():
            print('Error: Division failed.')
            assert(_sage_const_0 )
        if sB != R(_sage_const_0 ):
            print('Error: Remainder is not zero.')
            assert(_sage_const_0 )
        if gB.degree() > self.K_B.order or hB.degree() > max(pB.degree(), self.K_B.order - _sage_const_1  + qB.degree()):
            print('Error: Degree of gB or hB exceeds maximum bound.')
            assert(_sage_const_0 )

        ## C
        pC = self.H.vanishing_polynomial(x=gamma)*self.H.vanishing_polynomial(x=beta)*self.C.val()
        qC = (gamma - self.C.row())*(beta - self.C.col())
        points_C = []
        for k in self.K_C.to_list:
            points_C.append((F(k), (pC/qC)(x=k)))
        xgC = R.lagrange_polynomial(points_C) - omega_C / self.K_C.order
        gC, rC = xgC.quo_rem(R.lagrange_polynomial([(_sage_const_1 , _sage_const_1 ), (-_sage_const_1 , -_sage_const_1 )]))
        fC = xgC + omega_C / self.K_C.order
        if rC != R(_sage_const_0 ):
            print('Error: Remainder is not zero.')
            assert(_sage_const_0 )
        hC, sC = (pC - qC*fC).quo_rem(self.K_C.vanishing_polynomial())
        if pC - qC*fC != hC*self.K_C.vanishing_polynomial():
            print('Error: Division failed.')
            assert(_sage_const_0 )
        if sC != R(_sage_const_0 ):
            print('Error: Remainder is not zero.')
            assert(_sage_const_0 )
        if gC.degree() > self.K_C.order or hC.degree() > max(pC.degree(), self.K_C.order - _sage_const_1  + qC.degree()):
            print('Error: Degree of gC or hC exceeds maximum bound.')
            assert(_sage_const_0 )

        output_elements_to_file['hA'] = hA
        output_elements_to_file['hB'] = hB
        output_elements_to_file['hC'] = hC
        output_elements_to_file['gA'] = gA
        output_elements_to_file['gB'] = gB
        output_elements_to_file['gC'] = gC

        return (hA, hB, hC, gA, gB, gC)

    # NOTE: I changed this from our spec.
    def Round_6_lhs(self, hs, deltas):
        hA = hs[_sage_const_0 ]
        hB = hs[_sage_const_1 ]
        hC = hs[_sage_const_2 ]

        delta_A = deltas[_sage_const_0 ]
        delta_B = deltas[_sage_const_1 ]
        delta_C = deltas[_sage_const_2 ]

        h2 = delta_A * self.K_A.selector * hA * self.K_A.vanishing_polynomial()
        h2 += delta_B * self.K_B.selector * hB * self.K_B.vanishing_polynomial()
        h2 += delta_C * self.K_C.selector * hC * self.K_C.vanishing_polynomial()

        h2, r2 = h2.quo_rem(self.K.vanishing_polynomial()) # divide through by v_K

        output_elements_to_file['h2'] = h2

        return h2