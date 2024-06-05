using Oscar
using Combinatorics
pm = Polymake

currentDir = pwd()
include(joinpath(currentDir, "src/inputData38.jl"));
include(joinpath(currentDir, "src/fileHandling.jl"));
include(joinpath(currentDir, "src/tscCoordRing.jl"));
include(joinpath(currentDir, "src/matroidalSubd.jl"));
include(joinpath(currentDir, "src/Bmaximal.jl"));
include(joinpath(currentDir, "src/simplifyIdeal.jl"));
include(joinpath(currentDir, "G5Data/finDecompositionFunctions"));

G5Path = joinpath(currentDir,"groupsFinal/G5.dat")
G5 = file2SetVectors(G5Path);

R, x =  makePolyRing(3, 8, QQ)
dims_G5 = []

dimFile = joinpath(currentDir, "G5Data/dims_G5.dat")
io1 = open(dimFile, "w") 
close(io1)


for i in 1:length(G5)    
    w = G5[i]
    d = dimLimitFinCVP(w, QQ, 3, 8, R, x)
    push!(dims_G5, d)
    
    open(dimFile, "a") do io
        write(io, string("dim w", i, " = " , d, "\n"))
        end
    
end
