from sage.all import *
# Generates R1CS instances of m x n matrices where z is of the form [1, b^(d+2), b, b^2, ..., b, b, ...] and d is the multiplicative depth of the circuit
# this follows TestCircuit::gen_rand(..)
def gen_r1cs_instance(m, n, b, d):
    # generate the (public and private) witness vector belonging to our TestCircuit
    # the private part
    dummy_vars = n - 3 - d
    w = [2, 4]
    for i in range(0, dummy_vars):
        w.append(w[0])
    # the public part
    x = [1, 8]
    for d in range(3, d+2):
        x.append(w[1] * x[-1])
    z = x + w

    z = vector(z)
    # Generate constraints as per TestCircuit::generate_constraints(...)
    # We initialize the matrix so we can append to it, and cut off the first row later using submatrix
    A = matrix(zero_vector(m))
    B = matrix(zero_vector(m))
    C = matrix(zero_vector(m))

    num_mul_constraints = d - 1 - 1
    # Insert constraints of the form z[1]*z[2]=z[0]
    # TODO: this is hardcoded, make it dynamic based on TestCircuit
    for i in range(0, n - num_mul_constraints):
        a = zero_vector(m)
        a[4] = 1
        A = A.insert_row(A.nrows(), a)

        b = zero_vector(m)
        b[5] = 1
        B = B.insert_row(B.nrows(), b)

        c = zero_vector(m)
        c[1] = 1
        C = C.insert_row(C.nrows(), c)

    # Insert constraints of the form z[i-1]*z[2]=z[i]
    z_index = 2
    for i in range(0, num_mul_constraints):
        a = zero_vector(m)
        a[z_index - 1] = 1
        A = A.insert_row(A.nrows(), a)

        b = zero_vector(m)
        b[5] = 1
        B = B.insert_row(B.nrows(), b)

        c = zero_vector(m)
        c[z_index] = 1
        C = C.insert_row(C.nrows(), c)

        z_index += 1

    # Take submatrices using submatrix(i,j,nr,nc), start at entry (i,j), use nr rows, nc cols
    A = A.submatrix(1, 0, n, m)
    B = B.submatrix(1, 0, n, m)
    C = C.submatrix(1, 0, n, m)
    if (A*z).pairwise_product(B*z) != C*z:
        print('Error: Invalid R1CS instance.')
        assert(0)

    return (A, B, C, z, w, x)
