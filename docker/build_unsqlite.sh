#!/bin/sh


echo "Building unqlite"

if [ -f "/usr/include/unqlite.h" ]; then
    exit 0;
fi


cd ~
mkdir build_unql
cd build_unql
git clone https://github.com/symisc/unqlite.git
mkdir build
cmake unqlite -B build -DCMAKE_INSTALL_PREFIX=/usr
cd build
make
make install
