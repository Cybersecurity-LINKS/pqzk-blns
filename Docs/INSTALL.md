# Requirements
The following dependencies must be installed:
- [OpenSSL](https://www.openssl.org/source/) 3.3.0
- [GMP](https://gmplib.org/) 6.2.0
- [NTL](https://libntl.org/) 11.5.1

NOTE: It may work with greater (or lesser) versions.
Tested with Ubuntu 22.04.

## OpenSSL
For simplicity, the following commands download and install OpenSSL in ```$HOME/openssl-3.3.0``` folder:
```sh
wget https://www.openssl.org/source/openssl-3.3.0.tar.gz
tar xf openssl-3.3.0.tar.gz
cd openssl-3.3.0
./Configure --release --prefix=$HOME/openssl-3.3.0 -Wl,--enable-new-dtags,-rpath,'$(LIBRPATH)'
make -j10
make test -j10
make install -j10

cd ..
```

Once OpenSSL is correctly installed, the command
```~/openssl-3.3.0/bin/openssl version```
should return:
```sh
$ OpenSSL 3.3.0 9 Apr 2024 (Library: OpenSSL 3.3.0 9 Apr 2024)
```

## GMP
Install GMP 6.2.0, since NTL does not work with GMP 6.3.0 (latest version). 
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
OpenSSL
- https://www.openssl.org/source/
- https://github.com/openssl/openssl/blob/master/INSTALL.md

GMP
- https://gmplib.org/

NTL
- https://libntl.org/doc/tour-unix.html
- https://libntl.org/doc/tour-modules.html
