# Sage tests

You need to pass three arguments `m n b d` where:
- m is the number of constraints
- n is the number of variables
- b is the "base" for z (i.e. z = [b, b^2, b^3, ..., b^n] )
- d is the multiplicative depth of the circuit

Example usage:
```
sage varuna.sage 7 7 2 3
```

This will create a circuit with `padded_public_variables`: `[1, 8, 32, 128]` and `private_variables`: `[2, 4, 2]`.

This will print test vectors to a file, which should be moved into `snarkVM/algorithms/src/snark/varuna/resources`.