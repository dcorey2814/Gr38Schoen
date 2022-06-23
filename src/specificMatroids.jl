function matroidU(l1, l2, l3, l4)
    ls = [l1, l2, l3, l4]
    hyp = [vcat(l1,l2), vcat(l1,l3), vcat(l1,l4), vcat(l2,l3), vcat(l2,l4), vcat(l3,l4)]
    n_M = sum([length(l) for l in ls])
    return matroid_from_hyperplanes(hyp, n_M)
end

function matroidV(l1,l2,l3,l4,l5)
    ls = [l1,l2,l3,l4,l5]
    hyp = [vcat(l1,l2,l3), vcat(l4,l5), vcat(l1,l4), vcat(l1,l5), vcat(l2,l4), vcat(l2,l5), vcat(l3,l4), vcat(l3,l5)]
    n_M = sum([length(l) for l in ls])
    return matroid_from_hyperplanes(hyp, n_M)
end


function matroidW(l1,l2,l3,l4,l5)
    ls = [l1,l2,l3,l4,l5]
    hyp = [vcat(l1,l2,l3), vcat(l1,l4,l5), vcat(l2,l4), vcat(l2,l5), vcat(l3,l4), vcat(l3,l5)]
    n_M = sum([length(l) for l in ls])
    return matroid_from_hyperplanes(hyp, n_M)
end
