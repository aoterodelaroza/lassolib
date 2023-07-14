#! /bin/bash

gfortran -pg -g -fbacktrace -fexternal-blas -fallow-argument-mismatch -std=legacy -c dgqvec.f
gfortran -pg -g -fbacktrace -fexternal-blas -fallow-argument-mismatch -std=legacy -c dqrtv1.f
gfortran -pg -g -fbacktrace -fexternal-blas -fallow-argument-mismatch -std=legacy -c dqrqh.f
gfortran -pg -g -fbacktrace -fexternal-blas -fallow-argument-mismatch -std=legacy -c dqhqr.f
gfortran -pg -g -fbacktrace -fexternal-blas -fallow-argument-mismatch -std=legacy -c dqrot.f
gfortran -pg -g -fbacktrace -fexternal-blas -fallow-argument-mismatch -std=legacy -c dqrinc.f
gfortran -pg -g -fbacktrace -fexternal-blas -fallow-argument-mismatch -std=legacy -c dqrinr.f
gfortran -pg -g -fbacktrace -fexternal-blas -c lassofun.f90
gfortran -pg -g -fbacktrace -fexternal-blas -c main.f90
gfortran -pg -g -fbacktrace -fexternal-blas -o main main.o lassofun.o dqrinc.o dqrinr.o \
	 dqrtv1.o dqrot.o dqrqh.o dqhqr.o dgqvec.o -lblas -llapack

# gfortran -fexternal-blas -fallow-argument-mismatch -std=legacy -c dgqvec.f
# gfortran -fexternal-blas -fallow-argument-mismatch -std=legacy -c dqrtv1.f
# gfortran -fexternal-blas -fallow-argument-mismatch -std=legacy -c dqrqh.f
# gfortran -fexternal-blas -fallow-argument-mismatch -std=legacy -c dqhqr.f
# gfortran -fexternal-blas -fallow-argument-mismatch -std=legacy -c dqrot.f
# gfortran -fexternal-blas -fallow-argument-mismatch -std=legacy -c dqrinc.f
# gfortran -fexternal-blas -fallow-argument-mismatch -std=legacy -c dqrinr.f
# gfortran -fexternal-blas -c lassofun.f90
# gfortran -fexternal-blas -c main.f90
# gfortran -fexternal-blas -o main main.o lassofun.o dqrinc.o dqrinr.o \
# 	 dqrtv1.o dqrot.o dqrqh.o dqhqr.o dgqvec.o -lblas -llapack

