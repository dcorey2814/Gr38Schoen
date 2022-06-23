##############
##   Input Data  ###
##############

S8 = Array(pm.load_data("group38"));
S8 = pm.to_one_based_indexing(S8);
Delta38 = Polymake.polytope.hypersimplex(3,8);
vDelta38 = Matrix{Int64}(Delta38.VERTICES[1:56,1:9]);
lin = transpose(vDelta38);
rays = Matrix{Int64}(Polymake.load_data("GrRays38.data"));
MaxConesStr = readlines("ConesDrOfGr38.data");
MaxCones38_0index = [[Tuple([parse(Int64, s) for s in split(ss,"#")]) for ss in split(C," ")] for C in MaxConesStr];
MaxCones38 = pm.to_one_based_indexing(MaxCones38_0index);
