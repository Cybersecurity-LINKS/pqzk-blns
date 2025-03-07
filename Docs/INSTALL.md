# Requirements
The following dependencies must be installed:
- [GMP](https://gmplib.org/) 6.3.0
- [NTL](https://libntl.org/) 11.5.1

NOTE: it may work with greater (or lesser) versions.

Tested with Ubuntu 22.04, where both GMP and NTL can be simply installed as:
```sh
sudo apt install libgmp10 libgmp-dev libntl44 libntl-dev
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
