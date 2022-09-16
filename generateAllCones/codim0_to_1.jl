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
allIndex2Ray = Dict{Int64, Vector{Int64}}([i => AllRays_input[i] for i in 1:length(AllRays_input)])
allRays = Dict{Vector{Int64}, Int64}([AllRays_input[i] => i for i in 1:length(AllRays_input)])

maxCones38Index = [coneRep2Cone(S8, raysSec, allRays, C) for C in MaxCones38];

codim1Raw = Dict{Vector{Int64}, Matrix{Int64}}()

for k in 1:length(MaxCones38)
    inCone = coneRep2ConeNoLineality(S8,rays,MaxCones38[k])
    (inRays,inFacets) = cone2RaysFacets(inCone); 
    for F in inFacets
        raysF = Matrix{Int64}([inRays[i,j] for i ∈ F, j ∈ 1:56 ])
        w = sumRowsPrimitive(raysF); 
        codim1Raw[w] = raysF; 
    end
end

codim_1 = distinctOrbits(S8,collect(keys(codim1Raw)))
interior2Rayscodim1 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim1Raw[v] for v in codim_1 ])
codim1Index = [rayMatrix2ConeIndex(C, allRays) for C in values(interior2Rayscodim1)]

listVectors2File(codim_1, "allRepsByCodim/codim_1.dat")
listVectors2File(codim1Index, "allRepsByCodim/codim1Index.dat")
