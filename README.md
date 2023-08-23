# Sage tests

You need to pass three arguments `n m b d` where:
- n is the number of variables
- m is the number of constraints
- b is the "base" for z (i.e. z = [b, b^2, b^3, ..., b^n] )
- d is the multiplicative depth of the circuit

Example usage (TODO: in TestCircuit this corresponds to num_constraints:5 and num_variables:4):
```
sage varuna.sage 5 5 2 1
```
