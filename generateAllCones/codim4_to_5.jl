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

codim4Index = file2SetVectors("allRepsByCodim/codim4Index.dat")
interior2Rayscodim4 = Dict{Vector{Int64}, Matrix{Int64}}([coneIndex2Pair(CI, allRays, 56) for CI in codim4Index])

codim5Raw = coneRaysUpCodim(interior2Rayscodim4)

codim_5 = distinctOrbits(S8,collect(keys(codim5Raw)))
interior2Rayscodim5 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim5Raw[v] for v in codim_5 ])
codim5Index = [rayMatrix2ConeIndex(C, allRays) for C in values(interior2Rayscodim5)]

listVectors2File(codim_5, "allRepsByCodim/codim_5.dat")
listVectors2File(codim5Index, "allRepsByCodim/codim5Index.dat")
