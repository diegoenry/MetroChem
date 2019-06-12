#!/bin/bash

src_list="AngleForce BondForce Math Minimize ReadPDB WritePDB Main"
obj_list=""
for src in ${src_list} ; do
    gfortran -c ${src}.f90
    obj_list="${obj_list} ${src}.o"
done 

gfortran -o water_optimize ${obj_list}


