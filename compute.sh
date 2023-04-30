#!/bin/sh

cd Mesh
Freefem++ mesh_generate.edp #mesh generate

cd ../PDE
gcc PDE.c -O2 -o PDE
./PDE 2 ../mesh/mesh01.msh 

cd ../描画
./plot 2 ../mesh/mesh01.msh PDE 3