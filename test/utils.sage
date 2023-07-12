# -*- coding: utf-8 -*-
"""
Test utils for running test vectors with Varuna
"""

import os
from datetime import datetime

# Write the results of the Varuna sagemath implementation to a file
def write_test_results_to_file(A, B, C, z, circuit=""):

    # Get the current date and time
    now = datetime.now()

    # Format it as a string: YYYY-MM-DD_HH-MM-SS
    date_str = now.strftime("%Y-%m-%d")
    now_str = now.strftime("%H:%M:%S")

    # Create outputs directory if it doesn't exist
    if circuit == "":
        output_dir = f"outputs/{date_str}/{now_str}"
        os.makedirs(output_dir, exist_ok=True)
    else:
        output_dir = f"outputs/{date_str}/{now_str}/{circuit}"
        os.makedirs(output_dir, exist_ok=True)

    # Write r1cs instance to test file
    with open(f"{output_dir}/r1cs.txt", 'w') as f_r1cs:
        f_r1cs.write('A')
        f_r1cs.write('\n')
        for i in range(0, len(A)):
            for j in range(0, len(A[0])):
                f_r1cs.write(str(A[i][j]) + ', ')
            f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('B')
        f_r1cs.write('\n')
        for i in range(0, len(B)):
            for j in range(0, len(B[0])):
                f_r1cs.write(str(B[i][j]) + ', ')
            f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('C')
        f_r1cs.write('\n')
        for i in range(0, len(C)):
            for j in range(0, len(C[0])):
                f_r1cs.write(str(C[i][j]) + ', ')
            f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('\n')
        f_r1cs.write('z')
        f_r1cs.write('\n')
        for i in range(0, len(z)):
            f_r1cs.write(str(z[i]) + ', ')

    f_r1cs.close()

    # Write randomness and test elements to file
    with open(f"{output_dir}/randomness.txt", 'w') as f_r:
        for key in randomness_to_file:
            value = randomness_to_file[key]
            f_r.write(str(key) + ' ' + str(value))
            f_r.write('\n')
            f_r.write('\n')

    f_r.close()

    with open(f"{output_dir}/outputs.txt", 'w') as f_t:
        for key in output_elements_to_file:
            value = output_elements_to_file[key]
            f_t.write(str(key) + ' ')
            for coeff in R(value):
                f_t.write(str(coeff) + ',')
            f_t.write('\n')
            f_t.write('\n')

    f_t.close()

    with open(f"{output_dir}/groups.txt", 'w') as f_g:
        for key in group_elements_to_file:
            group = group_elements_to_file[key]
            f_g.write(str(key) + ' ')
            for g in group:
                f_g.write(str(F(g))+ ',')
            f_g.write('\n')
            f_g.write('\n')

    f_g.close()

# Load varuna test vectors into the config
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

# Verify that Varuna has produced outputs which match the test vectors
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



def load_instance(circuit_path: str):
    def parse_matrix_section(lines, start_line):
        matrix = []
        for line in lines[start_line:]:
            if line.strip() in ['A', 'B', 'C']:
                break  # Stop if we reach the next matrix identifier
            row = [int(x) for x in line.strip().split(',') if x.strip()]
            if row:
                matrix.append(row)
        return matrix

    # Read the file content
    with open(circuit_path, 'r') as file:
        lines = file.readlines()

    # Identify where each matrix starts
    A_start = lines.index('A\n') + 1 if 'A\n' in lines else None
    B_start = lines.index('B\n') + 1 if 'B\n' in lines else None
    C_start = lines.index('C\n') + 1 if 'C\n' in lines else None

    # Parse matrices
    A = parse_matrix_section(lines, A_start) if A_start is not None else None
    B = parse_matrix_section(lines, B_start) if B_start is not None else None
    C = parse_matrix_section(lines, C_start) if C_start is not None else None

    return A, B, C


def load_witness(witness_path: str):
    def parse_array(line):
        return [int(x) for x in line.strip().strip('[]').split(',')]

    # Read the file content
    with open(witness_path, 'r') as file:
        lines = file.readlines()

    # Parse the arrays
    x = parse_array(lines[0])
    w = parse_array(lines[1])

    return x, w



