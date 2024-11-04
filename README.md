# JoRGS: <ins>Jo</ins>int <ins>R</ins>otation <ins>G</ins>gate <ins>S</ins>ynthesis

## Introduction
`JoRGS` is a synthesis tool for quantum rotation gates implemented in C/C++. 
Given a circuit containing commuting rotation gates, `JoRGS` helps to synthesize it into the Clifford+T gate set with minimal T-count.

## Build
To build the binary file `JoRGS`, just use the following:
```commandline
make
```

## Execution
The circuit format being simulated is `OpenQASM` used by IBM's [Qiskit](https://github.com/Qiskit/qiskit), and the gate set supported in this simulator now contains Rotation-X (rx), Rotation-Y (ry), Rotation-Z (rz), Rotation-XX (rxx), Rotation-YY (ryy), Rotation-ZZ (rzz), Phase (p), and Controlled-phase (cp).
One can find example circuits in the [examples](/examples) folder. 

The help message states the details:

```commandline
$ ./JoRGS --help
Options:
  --help                produce help message
  --in arg              qasm file string for synthesis
  --out arg             qasm file string after synthesis
  --prec arg (=30)      precision in bits (default: 30)
  --cost arg (=1000)    T-count of applying an independent single-gate rotation. Set it to a large number to disable applying single-gate rotations (default: 1000)
  --same                use Fourier state transformation for the same-angle special case

```

For example, the following command synthesizes [example/vqe_layer.qasm](/example/vqe_layer.qasm), which is a layer in a VQE circuit.
```commandline
./JoRGS --in examples/vqe_layer.qasm --out out.qasm --prec 30
```
Then the synthesized circuit is produced in `out.qasm`.
In the output circuit, two ancilla quantum registers are used.
It is assumed that the "add" register is initialized as 0, and the "frs" register has been initialized as the Fourier state.
Moreover, the method in [[C. Gidney, 2018]](https://quantum-journal.org/papers/q-2018-06-18-74/) can be applied for canceling the Toffoli gates, which is used to calculate the T-count, but we keep the original circuit for clarity.


For the case that all rotation gates have the same rotation angles, the `--same` argument can be applied.    
For example, the following command synthesizes [example/qaoa_layer.qasm](/example/qaoa_layer.qasm), which is a layer in a QAOA circuit.
The `--cost 44` parameter specifies the T-count of applying each single-gate rotation with the [RUS method](https://arxiv.org/abs/1311.1074) used in the Fourier state transform.
```commandline
./JoRGS --in examples/qaoa_layer.qasm --out out.qasm --prec 30 --cost 44 --same
```
Note that we do not consider the T-count of the inverse Fourier state transform, as stated in the paper, but we keep it in the output circuit for clarity.

