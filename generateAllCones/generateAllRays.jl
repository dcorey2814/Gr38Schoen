using Oscar
using Combinatorics
pm = Polymake;

currentDir = pwd()
include(joinpath(currentDir, "src/inputData38.jl"));
include(joinpath(currentDir, "src/fileHandling.jl"));
include(joinpath(currentDir, "src/matroidalSubd.jl"));
include(joinpath(currentDir, "src_generateAllCones/generateAllCones.jl"));


raysSec = rays[1:12,:]

allRays = Dict{Vector{Int64}, Int64}()

global ind = 1
for ri in 1:12
    r = raysSec[ri,:]
    for g in S8
        gr = r[g]
        if gr in keys(allRays)
            continue
        else
            allRays[gr] = ind
            global ind += 1
        end
    end
end


allRaysInverse = Dict{Int64, Vector{Int64}}([allRays[v] => v for v in keys(allRays)] )
allRays_List = [allRaysInverse[i] for i in 1:length(keys(allRays))]

listVectors2File(allRays_List, "allRepsByCodim/allRays.dat")
