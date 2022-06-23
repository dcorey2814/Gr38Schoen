##########################
##   matroidal data of subdivision   ##
##########################

function plueckerList2Matroid(L, d, n)
    dn = collect(powerset([1:n;], d, d))
    return matroid_from_bases(dn[L], 1:n) 
end

function subd2Matroids(subd, d, n)
    return [plueckerList2Matroid(r, d, n) for r in maximal_cells(subd) ]     
end


function edgesDualGraph(subd)
    Gra = subd.pm_subdivision.POLYHEDRAL_COMPLEX.DUAL_GRAPH
    gdual = Graphs.Graph{Graphs.Undirected}(Gra.ADJACENCY)
    edges = collect(Graphs.edges(gdual))
    
    maxCells = maximal_cells(subd)
    
    dualGraphEdgesAsVertices = [intersect(maxCells[Graphs.src(e)], maxCells[ Graphs.src(Graphs.reverse(e)) ] ) for e in edges]
    
    return dualGraphEdgesAsVertices
end


function subd2CodimLeq1Matroids(subd, d, n)
    maxMatroids = [plueckerList2Matroid(Vector{Int64}(pm.@convert_to Vector r), d, n) for r in maximal_cells(subd) ]    
    codim1Matroids = [plueckerList2Matroid(Vector{Int64}(pm.@convert_to Vector r), d, n) for r in edgesDualGraph(subd) ]  
    return (maxMatroids, codim1Matroids) 
end

function subd2DualGraph(subd)
    Gra = subd.pm_subdivision.POLYHEDRAL_COMPLEX.DUAL_GRAPH    
    graphOscar = Graphs.Graph{Graphs.Undirected}(Gra.ADJACENCY)
    return graphOscar
end

function subdMatroidsNonLeaves(subd, d, n)
    Mats = subd2Matroids(subd, d, n)
    Gra = subd.pm_subdivision.POLYHEDRAL_COMPLEX.DUAL_GRAPH
    graphOscar = Graphs.Graph{Graphs.Undirected}(Gra.ADJACENCY)
    notLeaves = [v for v in 1:Graphs.nv(graphOscar) if length(Graphs.all_neighbors(graphOscar, v)) > 1 ]
    return Mats[notLeaves]
end

function commonBasis(Ms)
    return intersect([bases(M) for M in Ms]...)
end

function subdMatroidsLeaves(subd, d, n)
    Mats = subd2Matroids(subd, d, n)
    Gra = subd.pm_subdivision.POLYHEDRAL_COMPLEX.DUAL_GRAPH
    graphOscar = Graphs.Graph{Graphs.Undirected}(Gra.ADJACENCY)
    Leaves = [v for v in 1:Graphs.nv(graphOscar) if length(Graphs.all_neighbors(graphOscar, v)) == 1 ]
    return Mats[Leaves]
end

function leavesGraph(graphOscar)
    return [v for v in 1:Graphs.nv(graphOscar) if length(Graphs.all_neighbors(graphOscar, v)) == 1 ]
end

function leafEdgesDualGraphAsCells(subd)
    Gra = subd.pm_subdivision.POLYHEDRAL_COMPLEX.DUAL_GRAPH
    graphOscar = Graphs.Graph{Graphs.Undirected}(Gra.ADJACENCY)
    edges = collect(Graphs.edges(graphOscar))
    lG = leavesGraph(graphOscar)
    leafEdges = [e for e in edges if (Graphs.src(e) in lG || Graphs.src(Graphs.reverse(e)) in lG )]
    maxCells = maximal_cells(subd)
    dualGraphEdgesAsVertices = [intersect(maxCells[Graphs.src(e)], maxCells[ Graphs.src(Graphs.reverse(e)) ] ) for e in leafEdges]
    return dualGraphEdgesAsVertices
end


function subd2LeavesAndCenter(subd, d, n)
    leafMaxMatroids = subdMatroidsLeaves(subd, d, n)
    leafEdgeMaxMatroids = [plueckerList2Matroid(Vector{Int64}(pm.@convert_to Vector r), d, n) for r in leafEdgesDualGraphAsCells(subd)]
    remainingMaxMatroids = subdMatroidsNonLeaves(subd, d, n)
    
    return (leafMaxMatroids, leafEdgeMaxMatroids, remainingMaxMatroids) 
end
