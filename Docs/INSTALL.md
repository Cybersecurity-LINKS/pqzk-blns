# Requirements
The following dependencies must be installed:
- [GMP](https://gmplib.org/) 6.3.0
- [NTL](https://libntl.org/) 11.5.1
- [Falcon](https://falcon-sign.info/) 2021-11-01 (Optional)
- [Clang](https://clang.llvm.org/) 14.0.0 (Optional - necessary if Falcon is used)

NOTE: it may work with greater (or lesser) versions.

Tested with Ubuntu 22.04, where GMP, NTL, and Clang can be simply installed as:
```sh
sudo apt install libgmp10 libgmp-dev libntl44 libntl-dev clang
```
To build and install them manually, please follow the instructions below.

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
If ```USE_FALCON``` flag is set to ```1``` in the [Makefile](../Makefile) (default), 
Falcon reference implementation is automatically downloaded and used by ```make```.

Otherwise, it can be manually downloaded and built as a library as follows.
```sh
wget https://falcon-sign.info/Falcon-impl-20211101.zip
unzip Falcon-impl-20211101.zip
cd Falcon-impl-20211101
make -j$(nproc)
ar rcs libfalcon.a codec.o common.o falcon.o fft.o fpr.o keygen.o rng.o shake.o sign.o vrfy.o
cd ..
```

# Download & Build
```sh
wget https://github.com/Cybersecurity-LINKS/pqzk-blns/archive/refs/heads/main.zip -O ./BLNS.zip
unzip ./BLNS.zip
cd  pqzk-blns-main

clear && make clean

make -j$(nproc) && ./BLNS
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