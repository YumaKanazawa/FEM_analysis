#!/bin/sh

# cd Mesh
# Freefem++ mesh_generate.edp #mesh generate

cd PDE 
echo "What is p?"
read Lp

rm "advection_diffusion_$Lp.txt"
gcc PDE.c -O2 -o PDE
./PDE 2 ../mesh/mesh32.msh &./PDE 2 ../mesh/mesh48.msh &./PDE 2 ../mesh/mesh64.msh &./PDE 2 ../mesh/mesh80.msh
 
# for ((n=32;n<=80;n+=16)); do
#     ./PDE 2 ../mesh/mesh$n.msh 
#     # echo "$i"
# done
cd ../

# cd ../描画
# ./plot 2 ../mesh/mesh01.msh PDE 3