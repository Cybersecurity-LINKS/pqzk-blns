BLNS Framework for PQ Anonymous Verifiable Credentials
======================================================

This repository contains an implementation of the framework for Post-Quantum (PQ) Anonymous Verifiable Credentials defined by Bootle, Lyubashevsky, Nguyen, and Sorniotti (BLNS) and available in https://eprint.iacr.org/2023/560.pdf


## Requirements

- [GMP](https://gmplib.org/) 6.3.0
- [NTL](https://libntl.org/) 11.5.1
- [Falcon](https://falcon-sign.info/) 2021-11-01 (Optional)
- [Clang](https://clang.llvm.org/) 14.0.0 (Optional - necessary if Falcon is used)

NOTE: it is possible to set ```USE_FALCON = 1``` (default) in the [Makefile](./Makefile), 
to automatically download and use the ```Falcon_keygen``` and ```Falcon_GSampler``` from the [Falcon](https://falcon-sign.info/) reference implementation, for better performance.

Otherwise, with ```USE_FALCON = 0``` the ```NTRU_TrapGen``` and ```GSampler``` function defined in [Lattice.cc](./Lattice.cc) are used.

For installation instructions and additional information, see [/Docs/INSTALL.md](./Docs/INSTALL.md)

## Usage
After cloning this repository, build and run the ```BLNS``` executable as follows:
```sh
make -j$(nproc)

./BLNS
```

## Acknowledgements
Work done in collaboration with the Cryptography and Number Theory research group ([CrypTO](https://crypto.polito.it/)) at the Politecnico di Torino,
in the framework of the [QUBIP](https://qubip.eu/) project.
