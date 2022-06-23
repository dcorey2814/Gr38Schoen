#########################################################
##  checking B-maximality, and computing dimensions of thin Schubert cells.  ####
#########################################################

function dimBRing(M, xi)
    basesM = bases(M)
    VarB = [a for a in basesM if length(intersect(a,xi)) == 2  ] 
    return length(VarB)
end

function basis2DimB(M)
    dimBDict = Dict{Vector{Int64}, Int64}([b => dimBRing(M,b) for b in bases(M)])    
    return dimBDict
end

function basis2DimB_Mins(M)
    DM = basis2DimB(M)
    minB = min(values(DM)...)
    return Dict{Vector{Int64}, Int64}([b => DM[b] for b in keys(DM) if DM[b] == minB])  
end


function dimTSC(M, F, xi, R, x)
    stratumM = TSC(M, F, xi, R, x)
    dBM = dimBRing(M, xi)
    return dBM - codim(stratumM)
end


function dimTSC_Optimized(M, F, R, x)
    BDims = basis2DimB_Mins(M)
    xi = first(keys(BDims))
    
    stratumM = TSC(M, F, xi, R, x)
    dBM = dimBRing(M, xi)
    return dBM - codim(stratumM)
end

function isBMaximal(M, F, R, x)
    DM = basis2DimB_Mins(M)
    xi = first(keys(DM))
    minB = DM[xi]
    dimStratum = dimTSC(M,F,xi,R,x)

    return dimStratum == minB
end

function dimsTSC_Ms(Ms, F, R, x)
    return Dict([M => dimTSC_Optimized(M, F, R, x) for M in Ms])
end


function checkBMaximal(Ms, dimMs, F, xi)
    #dimsTSC_Ms = dimsTSC_Ms(Ms, F)
    dimsBRing_Ms = Dict([M => dimBRing(M, xi) for M in Ms])
    return all([dimMs[M] == dimsBRing_Ms[M] for M in Ms])
end

function existsBasisBMaximal_DimsPrecomputed(Ms, dimMs, F, Bs)
    return [xi for xi in Bs if checkBMaximal(Ms, dimMs, F, xi)]
end
    
function existsBasisBMaximal_Ms(Ms, F, R, x, Bs)
    dimMs = dimsTSC_Ms(Ms, F, R, x)
    return existsBasisBMaximal_DimsPrecomputed(Ms, dimMs, F, Bs)
end


function isFinBMaximal(Fin, Ms, F, R, x)
    MsInFin = Dict([i=>Ms[i] for i in Fin[1]])
    MsExposedFin = Dict([i=>Ms[i] for i in Fin[1] if i ∉ Fin[2]])
    Bs = commonBasis(values(MsInFin))
    return existsBasisBMaximal_Ms(values(MsExposedFin), F, R, x, Bs)    
end

function allFinsBMaximal(Fins, Ms, F, R, x)
    return all([length(isFinBMaximal(Fin, Ms, F, R, x)) > 0 for Fin in Fins])
end


function bases2DimBDifference(M)
    D = basis2DimB(M)
    m = minimum(values(D))
    
    return Dict([b => D[b] - m for b in keys(D)])
end

function optimalBasesForLimit(Ms)
    DMs = Dict([M => bases2DimBDifference(M) for M in Ms])
    basesOptions = commonBasis(Ms)
    dimsByBasis = Dict([ b=> sum([DMs[M][b] for M in Ms]) for b in basesOptions])
    
    minDimSum = minimum(values(dimsByBasis))
    
    minDimsByBasis = [b for b in keys(dimsByBasis) if  dimsByBasis[b] == minDimSum ]
    return minDimsByBasis
end

function dimBRing_MatList(Ms, xi)
    basesMs = []
    for M in Ms
        append!(basesMs, bases(M))
    end
    unique!(basesMs)
    VarB = [a for a in basesMs if length(intersect(a,xi)) == 2  ] 
    return length(VarB)
end




function dimLimitTSC(Ms, F, R, x)

    optB = optimalBasesForLimit(Ms)
    stratumMs = limitTSC(Ms, F, first(optB), R, x)
    gensSMs = minimal_generating_set(stratumMs)
  
    expM = exponentVecs(gensSMs, R, x)
    dBMs = dimBRing_MatList(Ms, first(optB))

    if searchUniUpperTriangular(expM)
        result = dBMs - length(gensSMs)
    else
        result = dBMs - codim(stratumMs)
    end
    
    return result
end


