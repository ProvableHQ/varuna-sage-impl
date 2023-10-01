# Varuna for Sagemath: A Second Implementation of Aleo's Zero-Knowledge Proof System

## Introduction
Welcome to the Sagemath implementation of the Varuna proof system. 

Varuna is an extension of the well-known Marlin proof system and has been implemented in production by the Aleo in Rust
to power the zero-knowledge proofs within the Aleo protocol. 

## Test Vectors
Test vectors verify the correctness of the implementation can be found in the `test` directory. Within this directory 
are folders that each contain:
1. `instance.input`: The R1CS instance the used in the proof
2. `witness.input`: A valid witness to satisfying the R1CS instance the prover is proving knowledge of
3. `challenge.input`: Verifier challenges provided to the prover (generated from Random oracles)
4. `domain`: A folder contain the various domains used by the PIOPs in the proof
5. `polynomials`: The intermediate PIOPs output by the prover at each step

These test vectors are provided to verify the correctness of this implementation and can be used to verify any other
concrete implementations of Varuna.

## Usage

Varuna can be run on a default example of a circuit representing `x^3 + x + 5 = 35` by invoking

```
sage run.sage
```

### Run Varuna on the Test Vectors
Varuna can be run on the existing test vectors by selecting any circuit in the `test` directory. The circuits are all in
named folders so you can select any circuit's test vectors you wish to run Varuna on.
```
sage varuna.sage circuit_0
```

If Varuna's output matches the test vectors the a message indicating the test vectors are valid will be printed.

Otherwise an error message will be printed indicating the test vectors are invalid.

### Run Varuna on a Custom Circuit
To run Varuna on a custom circuit - you need to pass four arguments `m n b d` where:
- m is the number of constraints
- n is the number of variables
- b is the "base" for z (i.e. z = [b, b^2, b^3, ..., b^n] )
- d is the multiplicative depth of the circuit

Example usage:
```
sage run.sage 7 7 2 3
```

The example above will create a circuit with `padded_public_variables`: `[1, 8, 32, 128]` and `private_variables`: `[2, 4, 2]`.

### Outputs

All outputs are written to the `outputs` directory (which will be created upon invoking `run.sage` if doesn't exist. 
If you would like to see the results of the invocation of Varuna, please check the `outputs` directory.