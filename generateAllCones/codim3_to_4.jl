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

codim3Index = file2SetVectors("allRepsByCodim/codim3Index.dat")
interior2Rayscodim3 = Dict{Vector{Int64}, Matrix{Int64}}([coneIndex2Pair(CI, allRays, 56) for CI in codim3Index])

codim4Raw = coneRaysUpCodim(interior2Rayscodim3)

codim_4 = distinctOrbits(S8,collect(keys(codim4Raw)))
interior2Rayscodim4 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim4Raw[v] for v in codim_4 ])
codim4Index = [rayMatrix2ConeIndex(C, allRays) for C in values(interior2Rayscodim4)]

listVectors2File(codim_4, "allRepsByCodim/codim_4.dat")
listVectors2File(codim4Index, "allRepsByCodim/codim4Index.dat")
