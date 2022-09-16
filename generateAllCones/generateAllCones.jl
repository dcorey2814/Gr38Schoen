function coneRep2Cone(G, rayOrbs, allRs, C)
    Cindex = Vector{Int64}()
    for i in 1:length(C)
        s = C[i]
        r = rayOrbs[s[1],G[s[2]]]
        push!(Cindex, allRs[r])
    end
    return sort!(Cindex)
end

function rayMatrix2ConeIndex(rayM, allRs)
    rs, cs = size(rayM)
    Cindex = [allRs[rayM[i,:]] for i in 1:rs]
    return sort!(Cindex)
end

function coneIndex2Pair(CI, allRs, nc)
    rayM = [allIndex2Ray[i][j] for i in CI, j in 1:nc]
    v = sumRowsPrimitive(rayM)
    return (v, rayM)
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
