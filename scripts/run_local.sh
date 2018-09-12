#!/usr/bin/env bash
fileName="laplace"

mkdir -p ../out

mpicc ../src/main.c -o ../out/${fileName}

mpirun -np 2 ../out/${fileName}