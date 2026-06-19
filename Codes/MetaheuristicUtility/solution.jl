mutable struct TransitSolution
    lineplan
    obj::Float64
    shortest
    directlink
    path_old
    lineinfo
end

TransitSolution(lineplan) = TransitSolution(lineplan, Inf, nothing, nothing, nothing, nothing)