import os
import ast

work_dir = os.getcwd();
load_attach_path(work_dir + "/algebra")
load_attach_path(work_dir + "/snark")
resources_folder = work_dir + "/test_vectors"
iop_polynomials = ["g_1", "g_a", "g_b", "g_c", "h_0", "h_1", "h_2", "w_lde", "z_lde"]
iop_domains = ["C", "K"]
verifier_challenges = ["alpha", "eta_A", "eta_B", "eta_C", "beta", "delta_A", "delta_B", "delta_C", "gamma"]

load('snark/indexer.sage')
load("snark/r1cs.sage")
load("snark/prover.sage")
load("snark/verifier.sage")

test_vectors = {}
group_elements_to_file = {}
matrix_elements_to_file = {}
output_elements_to_file = {}
randomness_to_file = {}

def test_cases(A, B, C, z, w=None, x=None):

    I = Indexer(matrix(A), matrix(B), matrix(C), vector(z), vector(w), vector(x))
    P = Prover(I.A, I.B, I.C, I.K, I.K_A, I.K_B, I.K_C, I.variable_domain, I.constraint_domain, I.X, I.z, I.x_poly, I.w_poly)

    row_oracles = [P.A.row, P.B.row, P.C.row]
    col_oracles = [P.A.col, P.B.col, P.C.col]
    val_oracles = [P.A.val, P.B.val, P.C.val]

    V = Verifier(row_oracles, col_oracles, val_oracles, I.K, I.K_A, I.K_B, I.K_C, I.variable_domain, I.constraint_domain, I.X, P.z_poly, P.x_poly, P.w_poly)

    # PIOP 1: Rowcheck
    h = P.Round_1_lhs()
    (gamma, eta_A, eta_B, eta_C) = V.Round_1_rhs()
    (sigA, sigB, sigC) = P.Round_2_lhs(gamma)
    etas = [eta_A, eta_B, eta_C]
    bit_0 = V.Round_2_rhs(sigA, sigB, sigC, h, gamma)
    print('Result of Rowcheck: ', bit_0)

    # PIOP 2: Univariate sumcheck
    (sigma, h1, g1, sigma_A, sigma_B, sigma_C) = P.Round_3_lhs(gamma, etas)
    sigmas = [sigma_A, sigma_B, sigma_C]
    beta = V.Round_3_rhs()
    (omega_A, omega_B, omega_C) = P.Round_4_lhs(gamma, beta)
    omegas = [omega_A, omega_B, omega_C]

    bit_1 = V.Round_4_rhs(sigma, h1, g1, omegas, etas, beta)
    print('Result of Univariate sumcheck: ', bit_1)

    #PIOP 3: Ratsumcheck
    (hA, hB, hC, gA, gB, gC) = P.Round_5_lhs(omegas, gamma, beta)
    hs = [hA, hB, hC]
    gs = [gA, gB, gC]
    (deltaA, deltaB, deltaC) = V.Round_5_rhs()
    deltas = [deltaA, deltaB, deltaC]
    h2 = P.Round_6_lhs(hs, deltas)
    bit_2 = V.Round_6_rhs(gs, gamma, beta, deltas, omegas, h2)
    print('Result of Rational sumcheck: ', bit_2)

def main():
    args = sys.argv[1:]

    if args[0] == 'custom':
        # custom test case
        A = matrix([[0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 1, 0, 0], [5, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0]])
        B = matrix([[0, 1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]])
        C = matrix([[0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]])
        x = vector([1, 3, 35])
        w = vector([9, 27, 30, 0])
        z = vector([1, 3, 35, 9, 27, 30, 0])
        test_cases(A, B, C, z, w, x)
    else:
        circuit = ""
        m = int(args[0])
        n = int(args[1])
        b = int(args[2])
        d = int(args[3])
        if len(args) == 5:
            circuit = str(args[4])
            ## Ensure the circuit path exists
            circuit_path = resources_folder + "/" + circuit
            if not os.path.exists(circuit_path):
                print(f"Circuit at {circuit_path} does not exist - aborting test")
                return
            load_test_vectors(circuit_path)

        (A, B, C, z, w, x) = gen_r1cs_instance(m, n, b, d)
        A = matrix(A)
        B = matrix(B)
        C = matrix(C)
        test_cases(A, B, C, z, w, x)
        if len(test_vectors) > 0:
            if verify_test_vectors():
                print(f"\nVaruna verified for {circuit} test vectors!")
            else:
                print(f"\nVaruna failed to verify for {circuit} test vectors!")

    # Write r1cs instance to test file
    with open('r1cs.txt', 'w') as f_r1cs:
        f_r1cs.write('A')
        f_r1cs.write('\n')
        for i in range(0, A.nrows()):
            for j in range(0, A.ncols()):
                f_r1cs.write(str(A[i, j]) + ', ')
            f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('B')
        f_r1cs.write('\n')
        for i in range(0, B.nrows()):
            for j in range(0, B.ncols()):
                f_r1cs.write(str(B[i, j]) + ', ')
            f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('C')
        f_r1cs.write('\n')
        for i in range(0, C.nrows()):
            for j in range(0, C.ncols()):
                f_r1cs.write(str(C[i, j]) + ', ')
            f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('z')
        f_r1cs.write('\n')
        for i in range(0, len(z)):
            f_r1cs.write(str(z[i]) + ', ')

    f_r1cs.close()

    # Write randomness and test elements to file
    with open('randomness.txt', 'w') as f_r:
        for key in randomness_to_file:
            value = randomness_to_file[key]
            f_r.write(str(key) + ' ' + str(value))
            f_r.write('\n')
            f_r.write('\n')

    f_r.close()

    with open('outputs.txt', 'w') as f_t:
        for key in output_elements_to_file:
            value = output_elements_to_file[key]
            f_t.write(str(key) + ' ')
            for coeff in R(value):
                f_t.write(str(coeff) + ',')
            f_t.write('\n')
            f_t.write('\n')

    f_t.close()

    with open('groups.txt', 'w') as f_g:
        for key in group_elements_to_file:
            group = group_elements_to_file[key]
            f_g.write(str(key) + ' ')
            for g in group:
                f_g.write(str(F(g))+ ',')
            f_g.write('\n')
            f_g.write('\n')

    f_g.close()

# Load test vectors into the config
def load_test_vectors(circuit_path: str):
    print(f"Running test vectors for circuit at {circuit_path}")

    # Load challenges
    with open(f"{circuit_path}/challenges.input", "r") as file:
        challenge_data = file.readlines()
        for x in zip(verifier_challenges, challenge_data):
            test_vectors[x[0]] = F(int(x[1].strip()))

    # Load polynomials and domains
    for item in (("domain", iop_domains), ("polynomials", iop_polynomials)):
        (vector_folder, vector_names) = item
        vectors = [(vector_name, f"{circuit_path}/{vector_folder}/{vector_name}.txt") for vector_name in vector_names]
        for vector in vectors:
            (vector_name, vector_path) = vector
            with open(vector_path, 'r') as f:
                test_vector = ast.literal_eval(f.read().strip())
                vector = [F(int(x)) for x in test_vector]
                if vector_folder == "polynomials":
                    test_vectors[vector_name] = R(vector)
                else:
                    test_vectors[vector_name] = vector

def verify_test_vectors():
    non_matching_polynomials = []
    non_matching_domains = []
    for poly in iop_polynomials:
        try:
            output_elements_to_file[poly] == test_vectors[poly]
        except:
            non_matching_polynomials.append(poly)

    for domain in iop_domains:
        try:
            group_elements_to_file[domain] == test_vectors[domain]
        except:
            non_matching_domains.append(domain)

    if len(non_matching_polynomials) > 0 or len(non_matching_domains) > 0:
        print("Test vectors do not match")
        print(f"Non Matching Polynomials: {non_matching_polynomials}")
        print(f"Non Matching Domains: {non_matching_domains}")
        return False
    return True

if __name__ == "__main__":
    main()
