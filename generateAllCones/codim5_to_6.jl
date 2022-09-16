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

codim5Index = file2SetVectors("allRepsByCodim/codim5Index.dat")
interior2Rayscodim5 = Dict{Vector{Int64}, Matrix{Int64}}([coneIndex2Pair(CI, allRays, 56) for CI in codim5Index])

codim6Raw = coneRaysUpCodim(interior2Rayscodim5)

codim_6 = distinctOrbits(S8,collect(keys(codim6Raw)))
interior2Rayscodim6 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim6Raw[v] for v in codim_6 ])
codim6Index = [rayMatrix2ConeIndex(C, allRays) for C in values(interior2Rayscodim6)]

listVectors2File(codim_6, "allRepsByCodim/codim_6.dat")
listVectors2File(codim6Index, "allRepsByCodim/codim6Index.dat")
