import os
import ast

# Setup the environment
work_dir = os.getcwd();
load_attach_path(work_dir + "/algebra")
load_attach_path(work_dir + "/snark")
load_attach_path(work_dir + "/test")
resources_folder = work_dir + "/test"

# Import Varuna logic
load('snark/indexer.sage')
load("snark/r1cs.sage")
load("snark/prover.sage")
load("snark/verifier.sage")
load("test/utils.sage")

# Specify the names of the intermediate polynomials, iop domains, and verifier challenges used in the Varuna proof
iop_polynomials = ["g_1", "g_a", "g_b", "g_c", "h_0", "h_1", "h_2", "w_lde", "z_lde"]
iop_domains = ["C", "K"]
verifier_challenges = ["alpha", "eta_A", "eta_B", "eta_C", "beta", "delta_A", "delta_B", "delta_C", "gamma"]

# Initialize the data containers used to contain proof outputs and the test vectors
group_elements_to_file = {}
matrix_elements_to_file = {}
output_elements_to_file = {}
randomness_to_file = {}
test_vectors = {}

# Execute a Varuna Snark proof for an R1CS instance
def run_varuna(A, B, C, z, w=None, x=None):

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
    ## TODO - Refactor this to use argparse so that named arguments can be passed in.
    cli_inputs = sys.argv

    if len(cli_inputs) == 1:
        # Run the Varuna on the circuit describing x^3 + x+ 5 = 35.
        A = matrix([[0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0], [0, 1, 0, 0, 1, 0, 0], [5, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0]])
        B = matrix([[0, 1, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]])
        C = matrix([[0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0], [0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0]])
        x = vector([1, 3, 35])
        w = vector([9, 27, 30, 0])
        z = vector([1, 3, 35, 9, 27, 30, 0])
        run_varuna(A, B, C, z, w, x)
        write_test_results_to_file(A, B, C, z)
    elif len(cli_inputs) == 2:
        # Run the Varuna on the test vectors specified by the user.
        # TODO - Refactor to support for reading arbitrary/multiple circuits from instance.input.

        # Get the name of the circuit Varuna should be tested on.
        circuit = cli_inputs[1]
        if cli_inputs[1] == "circuit_0":
            # Ensure the path to the circuit Varuna should be tested on exists path exists.
            circuit_path = f"{resources_folder}/{circuit}"
            if not os.path.exists(circuit_path):
                print(f"Circuit at {circuit_path} does not exist - aborting test")
                return

            # Load the test vectors
            load_test_vectors(circuit_path)

            # Generate the R1CS that Varuna will run on.
            (A, B, C, z, w, x) = gen_r1cs_instance(7, 7, 2, 3)
            A = matrix(A)
            B = matrix(B)
            C = matrix(C)

            # Run Varuna on the R1CS.
            run_varuna(A, B, C, z, w, x)

            # Verify that the intermediate outputs of Varuna match the expected test vectors.
            if len(test_vectors) > 0:
                if verify_test_vectors():
                    print(f"\nVaruna verified for {circuit} test vectors!")
                else:
                    print(f"\nVaruna failed to verify for {circuit} test vectors!")
        else:
            print(f"No test vectors exist for circuit {circuit} - aborting test")
            return

        # Write test results to file.
        write_test_results_to_file(A, B, C, z, circuit)
    elif len(cli_inputs) == 5:
        # Run Varuna on a randomly generated R1CS instance.
        args = cli_inputs[1:]
        m = int(args[0])
        n = int(args[1])
        b = int(args[2])
        d = int(args[3])

        # Generate the R1CS that Varuna will run on.
        (A, B, C, z, w, x) = gen_r1cs_instance(m, n, b, d)
        A = matrix(A)
        B = matrix(B)
        C = matrix(C)

        # Run Varuna on the R1CS.
        run_varuna(A, B, C, z, w, x)

        # Write test results to file.
        write_test_results_to_file(A, B, C, z)
    else:
        print("Invalid number of arguments passed to Varuna - aborting test")

if __name__ == "__main__":
    main()
