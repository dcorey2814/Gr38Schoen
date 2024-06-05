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

passes_test = []
fails_test = []

isBMax = joinpath(currentDir,"G5Data/isBMax.dat")
isNotBMax = joinpath(currentDir,"G5Data/isNotBMax.dat")

io1 = open(isBMax, "w") 
io2 = open(isNotBMax, "w")
close(io1)
close(io2)

R, x = makePolyRing(3,8, QQ)

for w in G5
    subd = subdivision_of_points(vDelta38[:,2:9], -w); 
    Ms = subd2Matroids(subd, 3, 8)
    DM = fin_decomp(w)
    if allFinsBMaximal(DM["fin"], Ms, QQ, R, x)
        push!(passes_test, w)
        open(isBMax, "a") do io
           write(io, vec2String(w), "\n")
        end;
    else
        push!(fails_test, w)
        open(isNotBMax, "a") do io
           write(io, vec2String(w), "\n")
        end;
    end    
end
