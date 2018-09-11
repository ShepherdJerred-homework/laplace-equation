#!/usr/bin/env bash
fileName="laplace"

mpicc src/main.c -o out/${fileName}

mpirun -np 2 out/${fileName}