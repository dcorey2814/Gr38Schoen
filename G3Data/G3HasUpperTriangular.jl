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

hasUpperTriangularG3 = []
doesntHaveUpperTriangularG3 = []

R, x = makePolyRing(3,8, QQ)

UTFile = joinpath(currentDir,"G3Data/hasUpperTriangularG3.dat")
NUTFile = joinpath(currentDir,"G3Data/doesnthaveUpperTriangularG3.dat")

io1 = open(UTFile, "w") 
io2 = open(NUTFile, "w")
close(io1)
close(io2)

hasUpperTriangularG3 = []
doesntHaveUpperTriangularG3 = []

R, x = makePolyRing(3,8, QQ)

UTFile = joinpath(currentDir,"G3Data/hasUpperTriangularG3.dat")
NUTFile = joinpath(currentDir,"G3Data/doesnthaveUpperTriangularG3.dat")

io1 = open(UTFile, "w") 
io2 = open(NUTFile, "w")
close(io1)
close(io2)


for w in G3
    subd = subdivision_of_points(vDelta38[:,2:9], -w)
    Ms = subdMatroidsNonLeaves(subd,3,8)
    
    optBs = optimalBasesForLimit(Ms)
    Lim = limitTSC(Ms, QQ, first(optBs), R, x)
    
    Rt = base_ring(Lim)
    xt = gens(base_ring(Lim))
    gLim = minimal_generating_set(Lim)
    
    if systemHasUniUpperTriangle(gLim, Rt, xt)
        push!(hasUpperTriangularG3, w)

        open(UTFile, "a") do io
            write(io, vec2String(w), "\n")
        end;

    else
        push!(doesntHaveUpperTriangularG3, w)
        open(NUTFile, "a") do io
            write(io, vec2String(w), "\n")
        end;
    end
    
end
