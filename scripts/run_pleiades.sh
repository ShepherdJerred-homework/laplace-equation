#!/usr/bin/env bash
fileName="laplace"

mpicc src/main.c -o out/${fileName}

mpirun -np 64 -hostfile ~/hostfile out/${fileName}