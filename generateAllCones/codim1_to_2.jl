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

codim1Index = file2SetVectors("allRepsByCodim/codim1Index.dat")
interior2Rayscodim1 = Dict{Vector{Int64}, Matrix{Int64}}([coneIndex2Pair(CI, allRays, 56) for CI in codim1Index])

codim2Raw = coneRaysUpCodim(interior2Rayscodim1)

codim_2 = distinctOrbits(S8,collect(keys(codim2Raw)))
interior2Rayscodim2 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim2Raw[v] for v in codim_2 ])
codim2Index = [rayMatrix2ConeIndex(C, allRays) for C in values(interior2Rayscodim2)]

listVectors2File(codim_2, "allRepsByCodim/codim_2.dat")
listVectors2File(codim2Index, "allRepsByCodim/codim2Index.dat")
