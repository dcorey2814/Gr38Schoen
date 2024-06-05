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


hasUpperTriangularG5 = []
doesntHaveUpperTriangularG5 = []

UTFile = joinpath(currentDir,"G5Data/hasUpperTriangularG5.dat")
NUTFile = joinpath(currentDir,"G5Data/doesnthaveUpperTriangularG5.dat")

io1 = open(UTFile, "w") 
io2 = open(NUTFile, "w")
close(io1)
close(io2)

R, x = makePolyRing(3,8, QQ)

for w in G5
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
    
    if systemHasUniUpperTriangle(gLim, Rt, xt)
        push!(hasUpperTriangularG5, w)
        open(UTFile, "a") do io
            write(io, vec2String(w), "\n")
        end;
    else
        push!(doesntHaveUpperTriangularG5, w)
        open(NUTFile, "a") do io
            write(io, vec2String(w), "\n")
        end;
    end 
end
