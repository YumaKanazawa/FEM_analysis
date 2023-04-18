using LinearAlgebra 
using CSV,DataFrames

cd("/Users/nakamura/Desktop/有限要素解析/FEM_ana/拡散方程式")
A=Matrix(CSV.read("file/heat_matrix.csv",DataFrame,header=false))   


for t in 0:100
    b=Matrix(CSV.read("file/init_matrix"*string(t)*".csv",DataFrame,header=false))
    
end
