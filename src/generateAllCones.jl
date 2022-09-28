using Oscar
using Combinatorics
pm = Polymake

S8 = Array(pm.load_data("group38"));
S8 = pm.to_one_based_indexing(S8);
Delta38 = Polymake.polytope.hypersimplex(3,8);
vDelta38 = Matrix{Int64}(Delta38.VERTICES[1:56,1:9]);
rays = Matrix{Int64}(Polymake.load_data("GrRays38.data"));
lin = transpose(vDelta38);
MaxConesStr = readlines("ConesDrOfGr38.data");
MaxCones38_0index = [[Tuple([parse(Int64, s) for s in split(ss,"#")]) for ss in split(C," ")] for C in MaxConesStr];
MaxCones38 = pm.to_one_based_indexing(MaxCones38_0index);

subd1 = SubdivisionOfPoints(vDelta38[:,2:9], -rays[1,:])
subd2 = SubdivisionOfPoints(vDelta38[:,2:9], -rays[2,:])
subd3 = SubdivisionOfPoints(vDelta38[:,2:9], -rays[3,:])
subd4 = SubdivisionOfPoints(vDelta38[:,2:9], -rays[4,:])
subd5 = SubdivisionOfPoints(vDelta38[:,2:9], -rays[5,:])
subd6 = SubdivisionOfPoints(vDelta38[:,2:9], -rays[6,:])
subd7 = SubdivisionOfPoints(vDelta38[:,2:9], -rays[7,:])
subd8 = SubdivisionOfPoints(vDelta38[:,2:9], -rays[8,:])
subd9 = SubdivisionOfPoints(vDelta38[:,2:9], -rays[9,:])
subd10 = SubdivisionOfPoints(vDelta38[:,2:9], -rays[10,:])
subd11 = SubdivisionOfPoints(vDelta38[:,2:9], -rays[11,:])
subd12 = SubdivisionOfPoints(vDelta38[:,2:9], -rays[12,:])
subds = [subd1,subd2,subd3,subd4,subd5,subd6,subd7,subd8,subd9,subd10,subd11,subd12]
rays_reduced = map(subd ->  Vector{Int64}(subd.pm_subdivision.MIN_WEIGHTS), subds)
rays_reduced = [rays_reduced[i][j] for i in 1:12, j in 1:56]


function plueckerList2Matroid(L, d, n)
    dn = collect(powerset([1:n;], d, d))
    return pm.matroid.Matroid(N_ELEMENTS=n, BASES = pm.to_zero_based_indexing(dn[L]))
end

function subd2Matroids(subd, d, n)
    return [plueckerList2Matroid(Vector{Int64}(pm.@convert_to Vector r), d, n) for r in maximal_cells(subd) ]    
end

function subd2MatroidsBasesIndex(subd, d, n)
    dn = collect(powerset([1:n;], d, d))
    basesList = [Vector{Int64}(pm.@convert_to Vector r) for r in maximal_cells(subd) ]
    return sort!(basesList)    
end

function GOrbitOneRay(G,r)
    nG = size(G)
    S=Set{Vector{Int64}}([r[g] for g in G])
    return S
end

function GOrbitsRays(G,Rs)
    S = Set{Vector{Int64}}()
    nRows = size(rays)[1]
    for i in 1:nRows
        Gr = GOrbitOneRay(G,Rs[i,:]) 
        if length(intersect(S,Gr))==0
            push!(S,Rs[i,:])
        end            
    end
    return S
end

function coneRep2InteriorVector(G,Rs,C)
    GRs = [ Rs[s[1],G[s[2]]] for s in C]
    return sum(GRs)
end

function coneRep2Cone(G,Rs,Ls,C)
    dRs = size(Rs);
    dLs = size(Ls);
    GRs = zeros(Int64, length(C), dRs[2]);
    for i in 1:length(C)
        s = C[i]
        GRs[i,:] = Rs[s[1],G[s[2]]]
    end
    return positive_hull(vcat(GRs, Ls, -Ls))
end

function coneRep2ConeNoLineality(G,Rs,C)
    dRs = size(Rs);
    GRs = zeros(Int64, length(C), dRs[2]);
    for i in 1:length(C)
        s = C[i]
        GRs[i,:] = Rs[s[1],G[s[2]]]
    end
    return positive_hull(GRs)
end


function facetAsTopMatroids_d3n8(raysC, f)
    raysfList = [raysC[i,:] for i in f]; 
    ncols = length(raysfList[1]); nrows = length(f); 
    raysf = [raysfList[i][j] for i in 1:nrows, j in 1:ncols]
    fCone = positive_hull(raysf)
    w = fCone.pm_cone.REL_INT_POINT
    subd = SubdivisionOfPoints(vDelta38[:,2:9], -w)
    return subd2MatroidsBasesIndex(subd, 3, 8)
end

function cone2RaysFacets(C)
    rs = Matrix{Rational{Int64}}(C.pm_cone.RAYS)
    nrays, dimAmbientSpace = size(rs)
    
    intRs = zeros(Int64, nrays, dimAmbientSpace)
    for i in 1:nrays
        intRs[i,:] = vector2Primitive(rs[i,:])
    end   
    fsRaw = C.pm_cone.RAYS_IN_FACETS
    nfacets, nrays2 = size(fsRaw) 
    fs = [[j for j in 1:nrays if fsRaw[i,j]] for i in 1:nfacets]
    return (intRs, fs)
end

function coneRep2RaysFacets(G,Rs,Crep)
    inC = coneRep2ConeNoLineality(G,Rs,Crep)
    return cone2RaysFacets(inC)
end
    
function sumRows(M)
    nrows,ncols = size(M)
    rM = [M[i,:] for i in 1:nrows ]
    return(sum(rM))
end

function sumRowsPrimitive(M)
    v = sumRows(M)
    denom = gcd(v)
    return Vector{Int64}(v/denom)
end

function vector2Primitive(v)
    denom = gcd(v)
    return Vector{Int64}(v/denom)
end

function raysFacets2InteriorPoints(Rs,Fs)
    return map(F -> sumRowsPrimitive([Rs[i,j] for i ∈ F, j ∈ 1:56 ]), Fs)
end

function vec2String(v)
    vecString = map(i -> string(i), v)
    return join(vecString, " ")
end


function listStrings2File(L, fileName)
    io = open(fileName,"w")
    newLines = join(L, "\n")
    write(io, newLines)
    close(io)
    return fileName
end

function listVectors2File(L, fileName)
    toListStrings = map(v -> join(map(i->string(i), v), " "), L)
    return listStrings2File(toListStrings, fileName)
end

function orbitVector(G,w)
    return Set{Vector{Int64}}([[w[g[i]] for i in 1:length(w)] for g in G])
end

function orbitSize2Reps(G,ws)
    size2Reps = Dict{Int64, Vector{Vector{Int64}}}()
    
    for w in ws
        orbw = orbitVector(G,w); 
        orbwSize = length(orbw); 

        if orbwSize ∉ keys(size2Reps)
            size2Reps[orbwSize] = [w]; 
            continue
        end
        if orbwSize ∈ keys(size2Reps)
            append!(size2Reps[orbwSize], [w]) 
        end
        
    end
    return size2Reps    
end


function multisetSupport(v)
    return sort(v)
end


function multisetSupp2Vector(ws)
    supp2Rep = Dict{Vector{Int64}, Vector{Vector{Int64}}}()
    
    for w in ws
        suppw = multisetSupport(w); 

        if suppw ∉ keys(supp2Rep)
            supp2Rep[suppw] = [w]; 
            continue
        end
        if suppw ∈ keys(supp2Rep)
            append!(supp2Rep[suppw], [w]) 
        end
        
    end
    return supp2Rep    
end

function distinctOrbits(G,ws)
    suppDict = multisetSupp2Vector(ws)
    reps = Vector{Vector{Int64}}()
    
    for sw in keys(suppDict)
        
        if length(suppDict[sw]) == 1
            append!(reps, suppDict[sw])
            continue
        end      
        
        rep2Orbit = Dict{Vector{Int64}, Set{Vector{Int64}}}()
        for w in suppDict[sw]             
            if all([w ∉ orb for orb in values(rep2Orbit) ])
                rep2Orbit[w] = orbitVector(G,w)
            end
            
        end
        append!(reps, collect(keys(rep2Orbit)))
    end
    return reps    
end


function orbitDict2File(D, fileName)
    ws = collect(keys(D))
    return listVectors2File(ws, fileName)
end

function coneRaysUpCodim(higherDim)
    interior2Rays = Dict{Vector{Int64}, Matrix{Int64}}()
    for w in keys(higherDim)
        wCone = positive_hull(higherDim[w])
        (wRays,wFacets) = cone2RaysFacets(wCone); 
        for F in wFacets
            raysF = Matrix{Int64}([wRays[i,j] for i ∈ F, j ∈ 1:56 ])
            v = sumRowsPrimitive(raysF); 
            interior2Rays[v] = raysF; 
        end
    end
    return interior2Rays
end

function makeOrbitFiles(G, reps, folderName)
    for i in 1:length(reps)
        orb_i = collect(orbitVector(S8, reps[i]))
        listVectors2File(orb_i, folderName*"/"*string(i)*".dat")
    end
    return folderName
end

function string2Int64Vector(s)
    return map(i -> parse(Int64, i), split(s))
end

function file2SetVectors(fileName)
    return map(s -> string2Int64Vector(s), readlines(fileName))
end


codim1Raw = Dict{Vector{Int64}, Matrix{Int64}}()
#for i in 1:1000
for i in 1:length(MaxCones38)
    inCone = coneRep2ConeNoLineality(S8,rays,MaxCones38[i])
    (inRays,inFacets) = cone2RaysFacets(inCone); 
    for F in inFacets
        raysF = Matrix{Int64}([inRays[i,j] for i ∈ F, j ∈ 1:56 ])
        w = sumRowsPrimitive(raysF); 
        codim1Raw[w] = raysF; 
    end
end


#codim_1 = distinctOrbits(S8,collect(keys(codim1Raw)))
#listVectors2File(codim_1, "codim_1.dat")

codim_1 = file2SetVectors("codim_1.dat")
interior2Rayscodim1 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim1Raw[v] for v in codim_1 ])

print("1\n")

codim2Raw = coneRaysUpCodim(interior2Rayscodim1)

#codim_2 = distinctOrbits(S8,collect(keys(codim2Raw)))
#listVectors2File(codim_2, "codim_2.dat")

codim_2 = file2SetVectors("codim_2.dat")
interior2Rayscodim2 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim2Raw[v] for v in codim_2 ])

print("2\n")

codim3Raw = coneRaysUpCodim(interior2Rayscodim2)

#codim_3 = distinctOrbits(S8,collect(keys(codim3Raw)))
#listVectors2File(codim_3, "codim_3.dat")

codim_3 = file2SetVectors("codim_3.dat")
interior2Rayscodim3 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim3Raw[v] for v in codim_3 ])

print("3\n")

codim4Raw = coneRaysUpCodim(interior2Rayscodim3)

#codim_4 = distinctOrbits(S8,collect(keys(codim4Raw)))
#listVectors2File(codim_4, "codim_4.dat")

codim_4 = file2SetVectors("codim_4.dat")
interior2Rayscodim4 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim4Raw[v] for v in codim_4 ])

print("4\n")

codim5Raw = coneRaysUpCodim(interior2Rayscodim4)

codim_5 = distinctOrbits(S8,collect(keys(codim5Raw)))
listVectors2File(codim_5, "codim_5.dat")

codim_5 = file2SetVectors("codim_5.dat")
interior2Rayscodim5 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim5Raw[v] for v in codim_5 ])

print("5\n")

codim6Raw = coneRaysUpCodim(interior2Rayscodim5)

codim_6 = distinctOrbits(S8,collect(keys(codim6Raw)))
listVectors2File(codim_6, "codim_6.dat")

codim_6 = file2SetVectors("codim_6.dat")
interior2Rayscodim6 = Dict{Vector{Int64}, Matrix{Int64}}([v=>codim6Raw[v] for v in codim_6 ])

print("6\n")

codim7Raw = coneRaysUpCodim(interior2Rayscodim6)

codim_7 = distinctOrbits(S8,collect(keys(codim7Raw)))
listVectors2File(codim_7, "codim_7.dat")




