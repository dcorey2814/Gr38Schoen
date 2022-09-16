using Oscar
using Combinatorics
pm = Polymake;

currentDir = pwd()
include(joinpath(currentDir, "src/inputData38.jl"));
include(joinpath(currentDir, "src/fileHandling.jl"));
include(joinpath(currentDir, "src/matroidalSubd.jl"));
include(joinpath(currentDir, "generateAllCones/generateAllCones.jl"));


raysSec = rays[1:12,:]


AllRays_input = file2SetVectors("allRepsByCodim/allRays.dat")

AllRaysMatrix = [AllRays_input[i][j] for i in 1:length(AllRays_input), j in 1:56 ]

allIndex2Ray = Dict{Int64, Vector{Int64}}([i => AllRays_input[i] for i in 1:length(AllRays_input)])
allRays = Dict{Vector{Int64}, Int64}([AllRays_input[i] => i for i in 1:length(AllRays_input)])

maxCones38Index = [coneRep2Cone(S8, raysSec, allRays, C) for C in MaxCones38];
maxConesRelativeInterior = [sumRowsPrimitive(AllRaysMatrix[MC,:]) for MC in maxCones38Index ]

listVectors2File(maxConesRelativeInterior, "allRepsByCodim/codim_0-test.dat")
listVectors2File(maxCones38Index, "allRepsByCodim/codim0Index-test.dat")
