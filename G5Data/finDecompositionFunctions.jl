
function not_cell(Mp,i)
    C = Set{Int64}()
    for j in 1:size(Mp,1)
        if j != i
            Mj = pm.row(Mp,j)
            C = union(C,Mj)
        end
    end
    return C
end

function get_edges(subd)
    dg = subdDualGraph(subd)
    edges_dg = collect(edges(dg))
    return Set([Set([src(e), dst(e)]) for e in edges_dg])
end

function find_fin(Mp, subd)
    edges_set = get_edges(subd)
#eliminate leaves
    not_leaves = Set{Int64}()
    for t in 1:size(Mp,1)
        Mt = pm.row(Mp,t)
        if length(Mt)>2
            push!(not_leaves,t)
        else
            int = intersect(Mt,not_cell(Mp,t))
            if length(int)>1
                push!(not_leaves,t)
            end
        end
    end
#collect fins and connecting edge
    fins_n_edges = Vector()
    for i in not_leaves
        body = Set{Int64}()
        for t in not_leaves
            if t != i
                Mt = pm.row(Mp,t)
                body = union(body,Mt)
            end
        end
        
        Mi = pm.row(Mp,i)
        int = intersect(Mi, body)
        
        if int in edges_set
            push!(fins_n_edges, (collect(Mi), collect(int)))
        end
    end 
    
    return(fins_n_edges)
end

function find_stable_edges(Mp, subd)
    finds = find_fin(Mp, subd)
    fins_with_stable_edges = Vector()
    edges_set = get_edges(subd)

    for i in 1:length(finds)
    
        fin = finds[i][1]
        finedge=finds[i][2]
        fin_edges = Set()
        stable_edges = Vector()
    
        for edge in edges_set
            nice_edge = collect(edge)
            if issubset(nice_edge,fin)
                push!(fin_edges,nice_edge)
            end
        end
    
        for edge in fin_edges
            if length(intersect(edge,finedge)) == 1
                push!(stable_edges,collect(edge))
            end
        end
        
        push!(fins_with_stable_edges,(fin,stable_edges))
    end
    
    return fins_with_stable_edges
end
#find_stable_edges(Mp,Gra,subd)

#find leaves
function find_leaves(Mp)
    #find leaves
    leaves = Vector()
    for t in 1:size(Mp,1)
        Mt = pm.row(Mp,t)
        if length(Mt)>2
        else
            int = intersect(Mt,not_cell(Mp,t))
            if length(int)==1
                push!(leaves,collect(Mt))
            end
        end
    end
    return(leaves)
end

#find body
function find_body(Mp, subd)
    finds = find_fin(Mp, subd)
    common_edges = Vector()
    for i in length(finds)
        push!(common_edges, finds[i][2])
    end
    
    edges_set =get_edges(subd)
    
#eliminate leaves
    not_leaves = Set{Int64}()
    for t in 1:size(Mp,1)
        Mt = pm.row(Mp,t)
        int = intersect(Mt,not_cell(Mp,t))
        if length(int)>1
            push!(not_leaves,t)
        end
    end
#remove fins
    not_fins = Vector{Int64}()
    for i in not_leaves
        body = Set{Int64}()
        for t in not_leaves
            if t != i
                Mt = pm.row(Mp,t)
                body = union(body,Mt)
            end
        end
        Mi = pm.row(Mp,i)
        int = intersect(Mi,body)
        if int in edges_set
        else
            push!(not_fins,i)
        end
    end 
    #collect vertices in definned component
    max_cells = Set()
    for edge in common_edges
        union!(max_cells, edge)
    end
    
    for i in not_fins
        Mi = pm.row(Mp, i)
        union!(max_cells, Mi)
    end
    max_cells_vec = collect(max_cells)
    return max_cells_vec
end
#find_body(Mp,subd)

function fin_decomp(w)
    subdcone = subdivision_of_points(vDelta38[:,2:9], -w)#subdivision corresponding to cone
    Gra = subdDualGraph(subd)
    Tmc = subdcone.pm_subdivision.TIGHT_SPAN
    Mp = Tmc.MAXIMAL_POLYTOPES

    edges_set = get_edges(subdcone)
    
    a=find_fin(Mp, subdcone)
    b=find_stable_edges(Mp, subdcone)
    c=find_leaves(Mp)
    d=find_body(Mp, subdcone)
    
    FD = Dict("fin" => a, "stable" => b, "leaves"=>c, "body" => d)
    return(FD)
end
    


function exposedEdgesFin(Fin, oG)
    exposedV = exposedVerticesFin(Fin)
    allV = Fin[1]
    N = Dict([v => [Set([v,i]) for i in all_neighbors(oG, v) if i in allV] for v in exposedV])
    edgesAsSets = union!(values(N)...) 
    return [sort!(collect(e)) for e in edgesAsSets]
end

function exposedVerticesFin(Fin)
    return [i for i in Fin[1] if i ∉ Fin[2]]    
end

function matroidOfEdge(Ms,e,n)
    B1 = bases(Ms[e[1]])
    B2 = bases(Ms[e[2]])
    return matroid_from_bases(intersect(B1,B2), n)
end;

function dimLimitFinCVP(w, F, d, n, R, x)
    D = fin_decomp(w)
    subd = subdivision_of_points(vDelta38[:,2:n+1], -w)
    oG = subdDualGraph(subd)
    Ms = subd2Matroids(subd,d,n)
        
    leafVMatroids = [Ms[i] for i in leavesGraph(oG)]
    leafEMatroids = [matroidOfEdge(Ms, e, n) for e in D["leaves"]]
    
    dimsLeavesV = Vector{Int64}([dimTSC_Optimized(M, F, R, x) for M in leafVMatroids]); 
    dimsLeavesE = Vector{Int64}([dimTSC_Optimized(M, F, R, x) for M in leafEMatroids]); 
    
    FinsExposedV = union!([exposedVerticesFin(Fin) for Fin in D["fin"]]...)
    FinsExposedE = union!([exposedEdgesFin(Fin, oG) for Fin in D["fin"]]...)
        
    FinsExposedVMats = [Ms[i] for i in FinsExposedV]
    FinsExposedEMats = [matroidOfEdge(Ms,e,n) for e in FinsExposedE]
        
    basesFins = [intersect([bases(Ms[i]) for i in  Fin[1]]...) for Fin in D["fin"]]
    
    FinsMats = [matroid_from_bases(Bs, n) for Bs in basesFins]
    
    dimsFinsV = Vector{Int64}([dimTSC_Optimized(M, F, R, x) for M in FinsExposedVMats]); 
    dimsFinsE = Vector{Int64}([dimTSC_Optimized(M, F, R, x) for M in FinsExposedEMats]); 
    dimsFinsF = Vector{Int64}([dimTSC_Optimized(M, F, R, x) for M in FinsMats ]);
    
    dimBody = dimLimitTSC([Ms[i] for i in D["body"]], F, R, x)
        
    return dimBody + sum(dimsLeavesV) - sum(dimsLeavesE) + sum(dimsFinsV) - sum(dimsFinsE) + sum(dimsFinsF)
    
end;
