using Oscar
using Combinatorics
pm = Polymake;

currentDir = pwd()
include(joinpath(currentDir, "src/inputData38.jl"));
include(joinpath(currentDir, "src/fileHandling.jl"));
include(joinpath(currentDir, "src/matroidalSubd.jl"));
include(joinpath(currentDir, "generateAllCones/generateAllCones.jl"));


AllRays_input = file2SetVectors("allRepsByCodim/allRays.dat")
allIndex2Ray = Dict{Int64, Vector{Int64}}([i => AllRays_input[i] for i in 1:length(AllRays_input)])
allRays = Dict{Vector{Int64}, Int64}([AllRays_input[i] => i for i in 1:length(AllRays_input)])

codim2Index = file2SetVectors("allRepsByCodim/codim2Index.dat")
interior2Rayscodim2 = Dict{Vector{Int64}, Matrix{Int64}}([coneIndex2Pair(CI, allRays, 56) for CI in codim2Index])

codim3Raw = coneRaysUpCodim(interior2Rayscodim2)

codim_3 = distinctOrbits(S8,collect(keys(codim3Raw)))
interior2Rayscodim3 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim3Raw[v] for v in codim_3 ])
codim3Index = [rayMatrix2ConeIndex(C, allRays) for C in values(interior2Rayscodim3)]

listVectors2File(codim_3, "allRepsByCodim/codim_3.dat")
listVectors2File(codim3Index, "allRepsByCodim/codim3Index.dat")
