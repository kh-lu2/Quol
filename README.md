# Quol 
Welcome to Quol - a smol quantum simulator, that uses gates and operators to create quantum circuits.

I wanted to make it as general as I could to be able to create circuits of differents shapes and sizes, so feel free to write your own quantum programs with it!

## How to Quol
### Matrix
Matrices are, well, matrices.

### Gate
Gates exist in context. Each Gate has a Matrix and specifies which qubits need to be transformed using it.

### Operator
Operators are several Gates. You can add an Operator to another Operator!

### Circuit
Circuits are generally lists of Gates, but you can add Operators to them, too.

## A simple quantum circuit
```cpp
Circuit circuit(2);
circuit.add_gate<HadamardMatrix>({0});
circuit.add_gate<CNOTMatrix>({0, 1});
```
Can you guess what it does?

## Bigger examples
### Shor
It's an implementation of Shor's algorithm for `n = 15, a = 7 and m = 16`. If you want to see more Shor, take a look at `classical_shor` branch, where you can find Shor's algorithm but computed purely classicaly!

### Grover
There are three example Grover's searches.

#### Grover 1
Performs search for x, such as `f(x) = 1`, where `f(x) = 1` for `x = 11`.

#### Grover 2
Performs search for x, such as `f(x) = 1`, where `f(x) = 1` for `x = 1, 5, or 8`.

#### Grover 3
Performs search for x, such as `f(x) = 11`, where `f(x)` swaps bits 1 and 3 of x and then nots bits 0 and 1.

