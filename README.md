# Varuna for SageMath: a Second Implementation of the Varuna Proof System

## Introduction
Welcome to the SageMath implementation of the Varuna proof system. 

Varuna is an optimized version of the well-known Marlin (https://eprint.iacr.org/2019/1047.pdf) proof system and currently powers zero-knowledge proofs within the Aleo Virtual Machine (snarkVM) protocol. 

A detailed protocol specification can be found in the `docs` directory. The protocol is implemented for a single circuit in R1CS with zero-knowledge. We leave the implementation of multi-circuit batching and lookups to future work. 

The purpose of this implementation is two-fold; first, to provide a high-level implementation of Varuna in SageMath for researchers and developers, and second, to sanity check the Varuna implementation in snarkVM written in Rust by validating it with test vectors. 

## Test Vectors
Test vectors can be found in the `test` directory. Within this directory are folders that represent R1CS circuits which 
contain:
1. `instance.input`: The R1CS instance the used in the proof.
2. `witness.input`: A valid witness to satisfying the R1CS instance the prover is proving knowledge of.
3. `challenge.input`: Verifier challenges provided to the prover.
4. `domain`: A folder that contains the various domains used by the PIOPs in the Varuna proof.
5. `polynomials`: The intermediate PIOPs output by the prover at each step.

Steps on running Varuna on the test vectors can be found in the [Usage](#usage) section.

## Usage
### Run Varuna on a Custom Circuit 
Varuna can be run on a custom circuit by selecting a pre-generated circuit in the `test` directory. An example usage is: 
```
sage run.sage circuit_0
```
If the Varuna protocol runs successfully, it will be indicated in the logs and the relevant outputs for the polynomial IOPs will be outputted in the directory `outputs`.

Otherwise, an error message will be logged indicating the test vectors are invalid.


### Outputs

All outputs are written to the `outputs` directory (which will be created upon invoking `run.sage`). 
If you would like to see the results of the invocation of Varuna, please check the `outputs` directory.
