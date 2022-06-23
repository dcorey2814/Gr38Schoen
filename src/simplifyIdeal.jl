#####################
###   Simplifying ideals  ####
#####################


function supportPolynomial(f, R, x)
    expVecs = [exponent_vector(f,i) for i in 1:length(f)]
    
    if length(expVecs) == 0
        return zeros(Int64,  length(x))
    end
    
    return sum(expVecs)
end

function supportSize(v)
    return length([a for a in v if a ≠ 0])
end

function isUnitVector(v)
    n = length(v)    
    Z = [i for i in 1:n if v[i] == Int64(0)]
    O = [i for i in 1:n if v[i] == Int64(1)]
    return (length(Z) == n-1 && length(O) == 1)
end


function addColumnStillUniUpperTriangular(nr, nc, pivotsM, optionsM)

    outp = []    
    notPivots = sort!([i for i in 1:nr if i ∉ pivotsM]); 

    
    for cindex in 1:length(optionsM)
        c = optionsM[cindex]
        cNonPivots = c[notPivots]
        if isUnitVector(cNonPivots)
            position1 = findall(x->x==1, cNonPivots)[1]
            
            newPivotsM = sort!(vcat(pivotsM, notPivots[position1]))
            
            newOptionsM = [optionsM[i] for i in 1:length(optionsM) if i ≠ cindex]
            push!(outp, (newPivotsM, newOptionsM) )
        else
            continue
        end
    end
    return outp    
end

function removeZeroRows(M)
    nr, nc = size(M)
    rowsM = [M[i,:] for i in 1:nr]
    nonZeroRows = [v for v in rowsM if supportSize(v) > 0]
    return Matrix{Int64}([nonZeroRows[i][j] for i in 1:length(nonZeroRows), j in 1:nc])
end


function searchUniUpperTriangular(M)
    Mno0 = removeZeroRows(M) # to define
    nr, nc = size(Mno0); #print("nr=", nr, " nc=", nc, "\n")
    
    if nr == 0
        return true
    end
    
    if nr == 1
        return any([a == Int64(1) for a in Mno0]) || all([a == Int64(0) for a in Mno0])
    end
    
    
    initialPivots = Vector{Int64}([])
    colsMno0 = [Mno0[:,i] for i in 1:nc]
    
    Stepa = addColumnStillUniUpperTriangular(nr, nc, initialPivots, colsMno0)
    #print("Stepa = ",Stepa, "\n")
    
    
    for i in 1:(nr-1)
        Stepb = []
        for i in 1:length(Stepa)
            addData = addColumnStillUniUpperTriangular(nr, nc, Stepa[i][1], Stepa[i][2])
            append!(Stepb, addData)
        end 
        Stepa = unique!(Stepb)

    end
    
    return length(Stepa) > 0
end


function exponentVecs(Igens, R, x)
    rowsEV = [supportPolynomial(f, R, x) for f in Igens]
    return Matrix{Int64}([rowsEV[i][j] for i in 1:length(rowsEV), j in 1:length(x)] )
end

function systemHasUniUpperTriangle(Igens, R, x)
    expVecs = exponentVecs(Igens, R, x)
    return searchUniUpperTriangular(expVecs)
end


function localizingSemiGroup(Ms, F, B, R, x)
    d = rank(Ms[1])
    n = length(matroid_groundset(Ms[1]))
    
    basesX = []

    for M in Ms        
        append!(basesX, projectedBases(M, F, B, R, x))
    end 
    
    sTotal = MPolyPowersOfElement(basesX[1])
    
    if length(basesX) == 1
        return sTotal
    end
    
    for i in 2:length(basesX)
        if basesX[i] ∉ sTotal
            sTotal = product(sTotal, MPolyPowersOfElement(basesX[i]))
        end
    end
    
    return sTotal
    
end


function localizeLimitRing(Ms, F, B, R, x)
    sTotal = localizingSemiGroup(Ms, F, B, R, x)
    return Localization(sTotal)    
end

function coefficientPoly(v, f, R)
    terms_v = [term(f,i) for i in 1:length(f) if v in vars(monomial(f,i))]
    
    
    if length(terms_v) == 0
        return R(0)
    else
        for t_v in terms_v
            dict_tv = Dict([pair for pair in factor(t_v)]) 
            if dict_tv[v] > 1
                return "can't isolate"
            else
                continue
            end
        end
    end
    
    fv = sum(terms_v)
    fvFactor = factor(fv)
    ufv = unit(fvFactor)
    fvFactorDict = Dict([pair for pair in fvFactor]) 
    coef_v = prod([k^(fvFactorDict[k]) for k in keys(fvFactorDict) if k ≠ v ])
    return ufv * coef_v
    
end

function varCanBeEliminated(v, f, R, S)
    cv = coefficientPoly(v, f, R)
    
    if typeof(cv) == String 
        if cv == "can't isolate"
            return false
        end
    end
    
    return cv in S
    
end


function solve4v(v, f, R, S, SR)
    
    if !varCanBeEliminated(v,f,R,S)
        return "can't solve"
    end
    
    cv = coefficientPoly(v, f, R)
    terms_notv = [term(f,i) for i in 1:length(f) if v ∉ vars(monomial(f,i))]
    
    if length(terms_notv) == 0
        return "can't solve"
    end
    
    return R(-1) * sum(terms_notv) // cv
end



function idealEliminateVariable(Igens, R, x, S, SR, yi, yj, fy)
    nr_x , nc_x = size(x); #print("nr_x = ", nr_x, "\n")
    nrnc_x = nr_x*nc_x; #print("nrnc_x = ", nrnc_x, "\n")
    Sx = [SR(x[i,j]) for i in 1:nr_x, j in 1:nc_x]; #print("Sx = ", Sx, "\n"); print("Sxij type = ", typeof(Sx[yi,yj]))
    Sx[yi,yj] = SR(fy); 
    
    #print("fy = ", fy, "\n")
    
    F = hom(R, SR, z->z, reshape(Sx, nrnc_x))
    
    newIgens = []
    for g in Igens
        #print("g = ", g, "\n")
        Fg = numerator(F(g))
        if Fg == 0
            continue
        end
        factoredFg = factor(Fg)
        factoredFgDict = Dict([pair for pair in factoredFg])
        newFactors = Dict()
        for h in keys(factoredFgDict)
            if h ∉ S
                newFactors[h] = factoredFgDict[h]
            end
        end
        new_g = prod([k^(newFactors[k]) for k in keys(newFactors)])
        push!(newIgens, new_g)        
    end
    
    J = ideal(R, newIgens)
    
    # update S
    
    return (unique!(gens(J)), R, x, S, SR)
end

function canReduceIdeal(Igens, R, x, S, SR)
    if systemHasUniUpperTriangle(Igens, R, x)
        return true
    end
    
    stop = false
    
    for f in Igens
        
        
        if stop
            continue
        end
        
        for v in vars(f)
            
            solve_v = solve4v(v, f, R, S, SR)
            
            if typeof(solve_v) == String
                if solve_v == "can't solve"
                    continue
                end
            else
                v_index = findall(t->t == v, x)[1]
                v_r = v_index[1]; #print("v_r = ", v_r, "\n")
                v_c = v_index[2]; #print("v_c = ", v_c, "\n")
                elim_v = idealEliminateVariable(Igens, R, x, S, SR, v_r, v_c, solve_v)
                stop = true
                
                return canReduceIdeal(elim_v[1], R, x, S, SR)
                
                
            end
        end      
    end
    
    if !stop
        return false
    end
    
end
