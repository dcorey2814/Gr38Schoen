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

doesntHaveUpperTriangularG3 = file2SetVectors(joinpath(currentDir, "G3Data/doesnthaveUpperTriangularG3-precomputed.dat"));

canFullySimplifyG3 = []
cannotFullySimplifyG3 = []

R, x = makePolyRing(3,8, QQ)

canFile = joinpath(currentDir,"G3Data/canFullySimplifyG3.dat")
cannotFile = joinpath(currentDir,"G3Data/cannotFullySimplifyG3.dat")

io1 = open(canFile, "w") 
io2 = open(cannotFile, "w")
close(io1)
close(io2)


for w in doesntHaveUpperTriangularG3
    subd = SubdivisionOfPoints(vDelta38[:,2:9], -w)
    Ms = subdMatroidsNonLeaves(subd,3,8)
    
    optBs = optimalBasesForLimit(Ms)
    Lim = limitTSC(Ms, QQ, first(optBs), R, x)
    
    Rt = base_ring(Lim)
    xt = gens(base_ring(Lim))
    gLim = minimal_generating_set(Lim)
    
    
    Sinv = localizingSemiGroup(Ms, QQ, first(optBs), R, x)
    SR = Localization(Sinv)[1]
    
    if canReduceIdeal(gLim, R, x, Sinv, SR)
        push!(canFullySimplifyG3, w)
        
        open(canFile, "a") do io
           write(io, vec2String(w), "\n")
        end;
    else
        push!(cannotFullySimplifyG3, w)
        
        open(cannotFile, "a") do io
           write(io, vec2String(w), "\n")
        end;
    end
end

