#!/bin/sh

rm -r CMakeFiles
rm cmake_install.cmake
rm Makefile
rm CMakeCache.txt

cmake ..
make 
