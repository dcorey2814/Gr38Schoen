using Oscar
using Combinatorics
pm = Polymake

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


G3Path = joinpath(currentDir,"groupsFinal/G3.dat")
G3 = file2SetVectors(G3Path);


function dimLimitLeafCVP(subd, F, d, n, R, x)
    (leafMaxMatroids, leafEdgeMaxMatroids, remainingMaxMatroids)  = subd2LeavesAndCenter(subd, d, n)
    
    dimsLeafMax = [dimTSC_Optimized(M, F, R, x) for M in leafMaxMatroids ]; 
    dimsLeafEdge = [dimTSC_Optimized(M, F, R, x) for M in leafEdgeMaxMatroids]
    dimCenter = dimLimitTSC(remainingMaxMatroids, F, R, x)
    
    return sum(dimsLeafMax) - sum(dimsLeafEdge) + dimCenter
end

R, x =  makePolyRing(3,8,QQ)
dims_G3 = []


dimFile = joinpath(currentDir, "G3Data/dims_G3.dat")
io1 = open(dimFile, "w") 
close(io1)


for i in 1:length(G3)
    w = G3[i]
    subd = SubdivisionOfPoints(vDelta38[:,2:9], -w)
    d = dimLimitLeafCVP(subd, QQ, 3, 8, R, x)
    push!(dims_G3, d)
    
    open(dimFile, "a") do io
        write(io, string("dim w", i, " = " , d, "\n"))
        end
    
end
