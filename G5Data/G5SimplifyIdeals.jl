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

canFullySimplifyG5 = []
cannotFullySimplifyG5 = []

R, x = makePolyRing(3,8, QQ)

canFullySimplifyG5 = []
cannotFullySimplifyG5 = []

R, x = makePolyRing(3,8, QQ)

canFile = joinpath(currentDir,"G5Data/canFullySimplifyG5.dat")
cannotFile = joinpath(currentDir,"G5Data/cannotFullySimplifyG5.dat")

io1 = open(canFile, "w") 
io2 = open(cannotFile, "w")
close(io1)
close(io2)

for w in doesntHaveUpperTriangularG5
    subd = subdivision_of_points(vDelta38[:,2:9], -w)
    FinDict = fin_decomp(w)
    Body = FinDict["body"]
    
    Ms = subd2Matroids(subd, 3, 8)
    MsBody = [Ms[i] for i in Body]
    
    optB = optimalBasesForLimit(MsBody)
    Lim = limitTSC(MsBody, QQ, first(optB), R, x)
    
    Rt = base_ring(Lim)
    xt = gens(base_ring(Lim))
    gLim = minimal_generating_set(Lim)
    
    Sinv = localizingSemiGroup(MsBody, QQ, first(optB), R, x)
    SR = localization(Sinv)[1]
    
    if canReduceIdeal(gLim, R, x, Sinv, SR)
        push!(canFullySimplifyG5, w)
        
        open(canFile, "a") do io
           write(io, vec2String(w), "\n")
        end;
    else
        push!(cannotFullySimplifyG5, w)
        
        open(cannotFile, "a") do io
           write(io, vec2String(w), "\n")
        end;
    end
end
