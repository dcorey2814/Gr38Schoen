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

function subdDualGraph(subd)
    Gra = subd.pm_subdivision.POLYHEDRAL_COMPLEX.DUAL_GRAPH
    E = pm.graph.edges(Gra)
    E = [[a for a in Set(e)] for e in Array(E)]
    E = pm.to_one_based_indexing(E)
    return graph_from_edges(E)
end

function edgesDualGraph(subd)
    dg = subdDualGraph(subd)
    dg_edges = collect(edges(dg))
    maxCells = maximal_cells(subd)
    dualGraphEdgesAsVertices = [intersect(maxCells[src(e)], maxCells[dst(e)]) for e in dg_edges]
    return dualGraphEdgesAsVertices
end


function subd2CodimLeq1Matroids(subd, d, n)
    maxMatroids = [plueckerList2Matroid(Vector{Int64}(pm.@convert_to Vector r), d, n) for r in maximal_cells(subd) ]    
    codim1Matroids = [plueckerList2Matroid(Vector{Int64}(pm.@convert_to Vector r), d, n) for r in edgesDualGraph(subd) ]  
    return (maxMatroids, codim1Matroids) 
end

function subdMatroidsNonLeaves(subd, d, n)
    Mats = subd2Matroids(subd, d, n)
    dg = subdDualGraph(subd)
    Leaves = leavesGraph(dg)
    notLeaves = [v for v in 1:n_vertices(dg) if !(v in Leaves)]
    return Mats[notLeaves]
end

function commonBasis(Ms)
    return intersect([bases(M) for M in Ms]...)
end

function subdMatroidsLeaves(subd, d, n)
    Mats = subd2Matroids(subd, d, n)
    dg = subdDualGraph(subd)
    Leaves = leavesGraph(dg)
    return Mats[Leaves]
end

function leavesGraph(dg)
    return [v for v in 1:n_vertices(dg) if length(all_neighbors(dg, v)) == 1 ]
end

function leafEdgesDualGraphAsCells(subd)
    dg = subdDualGraph(subd)
    edges_dg = collect(edges(dg))
    lG = leavesGraph(dg)
    leafEdges = [e for e in edges_dg if (src(e) in lG || dst(e) in lG )]
    maxCells = maximal_cells(subd)
    dualGraphEdgesAsVertices = [intersect(maxCells[src(e)], maxCells[dst(e)]) for e in leafEdges]
    return dualGraphEdgesAsVertices
end


function subd2LeavesAndCenter(subd, d, n)
    leafMaxMatroids = subdMatroidsLeaves(subd, d, n)
    leafEdgeMaxMatroids = [plueckerList2Matroid(Vector{Int64}(pm.@convert_to Vector r), d, n) for r in leafEdgesDualGraphAsCells(subd)]
    remainingMaxMatroids = subdMatroidsNonLeaves(subd, d, n)
    return (leafMaxMatroids, leafEdgeMaxMatroids, remainingMaxMatroids) 
end
