using Oscar
using Combinatorics
pm = Polymake

currentDir = pwd()
include(joinpath(currentDir, "src/inputData38.jl"));
include(joinpath(currentDir, "src/fileHandling.jl"));
include(joinpath(currentDir, "src/tscCoordRing.jl"));
include(joinpath(currentDir, "src/matroidalSubd.jl"));
include(joinpath(currentDir, "src/Bmaximal.jl"));
include(joinpath(currentDir, "src/simplifyIdeal.jl"));


G5Path = joinpath(currentDir,"groupsFinal/G5.dat")
G5 = file2SetVectors(G5Path);


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
    Gra = subd.pm_subdivision.POLYHEDRAL_COMPLEX.DUAL_GRAPH
    
    gdual = Graphs.Graph{Graphs.Undirected}(Gra.ADJACENCY)
    edges = collect(Graphs.edges(gdual))
    edges_set = Set{Set}()
    for i in 1:length(edges)
        Ei = edges[i]
        edge = Set([Ei.source, Ei.target])
        push!(edges_set,edge)
    end
    return(edges_set)
end

function find_fin(Mp,Gra,subd)
    
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
        int = intersect(Mi,body)
        
        if int in edges_set
            push!(fins_n_edges,(collect(Mi),collect(int)))
        end
    end 
    
    return(fins_n_edges)
end

function find_stable_edges(Mp,Gra,subd)
    finds = find_fin(Mp,Gra,subd)
    
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
function find_body(Mp,Gra,subd)
    
    finds = find_fin(Mp,Gra,subd)
   
    common_edges = Vector()
    
    for i in length(finds)
        push!(common_edges,finds[i][2])
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
            union!(max_cells,edge)
        end

    
        for i in not_fins
            Mi = pm.row(Mp,i)
            union!(max_cells, Mi)
        end
        max_cells_vec = collect(max_cells)
        return max_cells_vec
end
#find_body(Mp,subd)

function fin_decomp(w)
    subdcone = SubdivisionOfPoints(vDelta38[:,2:9], -w)#subdivision corresponding to cone

    Gra = subdcone.pm_subdivision.POLYHEDRAL_COMPLEX.DUAL_GRAPH
    Tmc = subdcone.pm_subdivision.TIGHT_SPAN
    Mp = Tmc.MAXIMAL_POLYTOPES
    Mc = subdcone.pm_subdivision.MAXIMAL_CELLS

    edges_set = get_edges(subdcone)
    
    a=find_fin(Mp,Gra,subdcone)
    b=find_stable_edges(Mp,Gra,subdcone)
    c=find_leaves(Mp)
    d=find_body(Mp,Gra,subdcone)
     
    FD = Dict("fin" => a, "stable" => b, "leaves"=>c, "body" => d)

return(FD)
end

canFullySimplifyG5 = []
cannotFullySimplifyG5 = []

R, x = makePolyRing(3,8, QQ)

canFile = joinpath(currentDir,"G5Data/canFullySimplifyG5.dat")
cannotFile = joinpath(currentDir,"G5Data/cannotFullySimplifyG5.dat")

io1 = open(canFile, "w") 
io2 = open(cannotFile, "w")
close(io1)
close(io2)


for w in doesntHaveUpperTriangularG5
    subd = SubdivisionOfPoints(vDelta38[:,2:9], -w)
    FinDict = fin_decomp(w)
    Body = FinDict["body"]
    
    Ms = subd2Matroids(subd, 3, 8)
    MsBody = [Ms[i] for i in Body]
    
    optB = optimalBasesForLimit(MsBody)
    Lim = limitTSC(MsBody, QQ, first(optB), R, x)
    
    Rt = base_ring(Lim)
    xt = gens(base_ring(Lim))
    gLim = minimal_generating_set(Lim)
    
    Sinv = localizingSemiGroup(MsBody, QQ, first(optB), R, x)
    SR = Localization(Sinv)
    
    if canReduceIdeal(gLim, R, x, Sinv, SR)
        push!(canFullySimplifyG5, w)
        
        open(canFile, "a") do io
           write(io, vec2String(w), "\n")
        end;
    else
        push!(cannotFullySimplifyG5, w)
        
        open(cannotFile, "a") do io
           write(io, vec2String(w), "\n")
        end;
    end
end

