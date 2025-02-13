# Requirements
The following dependencies must be installed:
- [GMP](https://gmplib.org/) 6.2.0
- [NTL](https://libntl.org/) 11.5.1

NOTE: It may work with greater (or lesser) versions.
Tested with Ubuntu 22.04.

## GMP
```sh
wget https://gmplib.org/download/gmp/gmp-6.2.0.tar.xz
tar xf gmp-6.2.0.tar.xz
cd gmp-6.2.0/
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

cd ..
```

# Download & Build
```sh
wget https://github.com/Cybersecurity-LINKS/blns/archive/refs/heads/main.zip -O ./BLNS.zip
unzip ./BLNS.zip
cd  blns-main

clear && make clean

make -j10 && ./BLNS
```

## Additional information
GMP
- https://gmplib.org/

NTL
- https://libntl.org/doc/tour-unix.html
- https://libntl.org/doc/tour-modules.html
