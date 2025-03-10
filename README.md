BLNS Framework for PQ Anonymous Verifiable Credentials
======================================================

This repository contains an implementation of the framework for Post-Quantum (PQ) Anonymous Verifiable Credentials defined by Bootle, Lyubashevsky, Nguyen, and Sorniotti (BLNS) and available in https://eprint.iacr.org/2023/560.pdf


## Requirements

- [GMP](https://gmplib.org/) 6.3.0
- [NTL](https://libntl.org/) 11.5.1

For installation instructions and additional information, see [./Docs/INSTALL.md](./Docs/INSTALL.md)

## Usage
After cloning this repository, build and run the ```BLNS``` executable as follows:
```sh
make -j$(nproc)

./BLNS
```

## Acknowledgements
Work done in collaboration with the Cryptography and Number Theory research group ([CrypTO](https://crypto.polito.it/)) at the Politecnico di Torino,
in the framework of the [QUBIP](https://qubip.eu/) project.
