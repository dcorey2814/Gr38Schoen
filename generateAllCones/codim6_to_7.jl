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

codim6Index = file2SetVectors("allRepsByCodim/codim6Index.dat")
interior2Rayscodim6 = Dict{Vector{Int64}, Matrix{Int64}}([coneIndex2Pair(CI, allRays, 56) for CI in codim6Index])

codim7Raw = coneRaysUpCodim(interior2Rayscodim6)

codim_7 = distinctOrbits(S8,collect(keys(codim7Raw)))
interior2Rayscodim7 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim7Raw[v] for v in codim_7 ])
codim7Index = [rayMatrix2ConeIndex(C, allRays) for C in values(interior2Rayscodim7)]

listVectors2File(codim_7, "allRepsByCodim/codim_7.dat")
listVectors2File(codim7Index, "allRepsByCodim/codim7Index.dat")
