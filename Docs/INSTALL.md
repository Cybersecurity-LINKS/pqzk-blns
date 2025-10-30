# Requirements
The following dependencies must be installed:
- [GMP](https://gmplib.org/) 6.3.0
- [NTL](https://libntl.org/) 11.5.1
- [Falcon](https://falcon-sign.info/) 2021-11-01 (Optional)
- [Clang](https://clang.llvm.org/) 18.1.3 (Optional - necessary if Falcon is used)

NOTE: it may work with greater (or lesser) versions.

Tested with Ubuntu 22.04 and 24.04, where the dependencies can be simply installed as:
```sh
sudo apt install wget tar xz-utils unzip make m4 g++ libgmp10 libgmp-dev libntl44 libntl-dev clang
```
Otherwise, please follow the instructions below to build and install them manually.

## GMP
```sh
wget https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz
tar xf gmp-6.3.0.tar.xz
cd gmp-6.3.0/
./configure
make
make check
sudo make install

cd ..
```

## NTL
```sh
wget https://libntl.org/ntl-11.5.1.tar.gz
gunzip ntl-11.5.1.tar.gz 
tar xf ntl-11.5.1.tar 
cd ntl-11.5.1/src/
./configure
make
make check
sudo make install

cd ../..
```

## Falcon (Optional)
The [Makefile](../Makefile) is configured with ```USE_FALCON = 1``` (default), 
to automatically download and use the [Falcon](https://falcon-sign.info/) reference implementation, for better performance.

It can also be manually downloaded and built as a library with the following commands:
```sh
wget https://falcon-sign.info/Falcon-impl-20211101.zip
unzip Falcon-impl-20211101.zip
cd Falcon-impl-20211101
make -j$(nproc)
ar rcs libfalcon.a codec.o common.o falcon.o fft.o fpr.o keygen.o rng.o shake.o sign.o vrfy.o
cd ..
```

# Download & Build
After installing the dependencies, you can download this repository, build and run the ```BLNS``` executable as follows:

```sh
wget https://github.com/Cybersecurity-LINKS/pqzk-blns/archive/refs/heads/main.zip -O ./BLNS.zip
unzip ./BLNS.zip
cd  pqzk-blns-main

make -j$(nproc) && ./BLNS
```

# Build & Run with Docker
As an alternative, it is possible to build and run the code with [Docker](https://docs.docker.com/), 
using [dockerfile_blns](../dockerfile_blns) that is located in the top-level directory of the repository.

First, build the Docker image from [dockerfile_blns](../dockerfile_blns) (it may take some time):
```sh
docker build -t blns_test -f .\dockerfile_blns .
```

Then, run the Docker container and the ```BLNS``` executable:
```sh
docker run -t blns_test ./BLNS
```

or open an interactive terminal to run the ```BLNS``` executable:
```sh
docker run -it blns_test /bin/bash

./BLNS
```


## Additional information
GMP
- https://gmplib.org/

NTL
- https://libntl.org/doc/tour-unix.html
- https://libntl.org/doc/tour-modules.html

Falcon
- https://falcon-sign.info/
- https://falcon-sign.info/falcon.pdf
- https://falcon-sign.info/impl/falcon.h.html

Docker 
- https://docs.docker.com/