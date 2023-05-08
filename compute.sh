#!/bin/sh

# cd Mesh
# Freefem++ mesh_generate.edp #mesh generate

cd PDE 
gcc PDE.c -O2 -o PDE
for ((i=1;i<=10;i++)); do
    ./PDE 2 ../mesh/mesh0$i.msh 
    # echo "$i"
done


# cd ../描画
# ./plot 2 ../mesh/mesh01.msh PDE 3