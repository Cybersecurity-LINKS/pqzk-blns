BLNS Framework for PQ Anonymous Verifiable Credentials
======================================================

This repository contains an implementation of the framework for Post-Quantum (PQ) Anonymous Verifiable Credentials defined by Bootle, Lyubashevsky, Nguyen, and Sorniotti (BLNS) and available in https://eprint.iacr.org/2023/560.pdf


## Requirements

- [GMP](https://gmplib.org/) 6.3.0
- [NTL](https://libntl.org/) 11.5.1
- [Falcon](https://falcon-sign.info/) 2021-11-01 (Optional)
- [Clang](https://clang.llvm.org/) 18.1.3 (Optional - necessary if Falcon is used)

NOTE: the [Makefile](./Makefile) is configured with ```USE_FALCON = 1``` (default), 
to automatically download and use the ```Falcon_keygen``` and ```Falcon_GSampler``` from the [Falcon](https://falcon-sign.info/) reference implementation, for better performance.

Otherwise, it is possible to set ```USE_FALCON = 0``` to use the ```NTRU_TrapGen``` and ```GSampler``` functions defined in [Lattice.cc](./Lattice.cc).

For installation instructions and additional information, see [/Docs/INSTALL.md](./Docs/INSTALL.md).

## Usage
After cloning this repository, build the ```BLNS``` and ```bench``` executables as follows:
```sh
make -j$(nproc)
```
then, run a demo of the ```BLNS``` framework:
```sh
./BLNS
```
and/or run the benchmarks ```bench```:
```sh
./bench
```

As an alternative, it is possible to build and run the code with Docker, using [dockerfile_blns](./dockerfile_blns).
For more details, see [/Docs/INSTALL.md](./Docs/INSTALL.md).

Note that parameters for different security levels (i.e., λ = 128, 192 or 256 bit) are available in the [params.h](./params.h) file.
The default configuration is ```lambda0``` = 128 bit, while to test 192 or 256 bit security levels you can uncomment the appropriate ```#define lambda0``` line at the beginning of [params.h](./params.h)

Moreover, the demo and benchmarks are configured to run with 16 attributes (see ```l0``` = 16 in [params.h](./params.h)),
6 of which are disclosed/revealed (see the variable ```idx_pub``` in [BLNS.cc](./BLNS.cc) and [bench.cc](./bench.cc)).


## Acknowledgements
Work done in collaboration with the Cryptography and Number Theory research group ([CrypTO](https://crypto.polito.it/)) at the Politecnico di Torino,
in the framework of the [QUBIP](https://qubip.eu/) project.
