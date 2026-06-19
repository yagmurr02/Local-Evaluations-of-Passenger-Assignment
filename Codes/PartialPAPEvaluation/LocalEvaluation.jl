include(joinpath(@__DIR__, "..", "CoreImplementation", "evaluator.jl"))
using DataStructures: PriorityQueue, dequeue!

# GROUND TRUTH RECONSTRUCTION
# Expand Floyd-Warshall "intermediate node" matrix into full stop sequence
function _expand_fw_stops(i::Int, j::Int, path::AbstractMatrix{Int})
    k = path[i, j]
    if k == 0
        return [i, j]
    else
        left  = _expand_fw_stops(i, k, path)
        right = _expand_fw_stops(k, j, path)
        return vcat(left[1:end-1], right)
    end
end

# Build (stops, lines_taken) from (path, lineinfo)
function _reconstruct_from_fw(o::Int, d::Int,
                              path::AbstractMatrix{Int},
                              lineinfo::AbstractMatrix{Int})
    stops = _expand_fw_stops(o, d, path)

    lines_taken = Int[]
    for t in 1:length(stops)-1
        u = stops[t]
        v = stops[t+1]
        push!(lines_taken, lineinfo[u, v])
    end

    return stops, lines_taken
end

# Ground truth affected OD pairs: cost changed OR stop sequence changed OR line sequence changed (even if stops same)
function affected_od_groundtruth(shortest_old::AbstractMatrix{Float64},
                                 shortest_new::AbstractMatrix{Float64},
                                 path_old::AbstractMatrix{Int},
                                 path_new::AbstractMatrix{Int},
                                 lineinfo_old::AbstractMatrix{Int},
                                 lineinfo_new::AbstractMatrix{Int};
                                 tol::Float64 = 1e-9)

    n = size(shortest_old, 1)
    A = Set{Tuple{Int,Int}}()

    @inbounds for o in 1:n, d in 1:n
        o == d && continue

        sold = shortest_old[o, d]
        snew = shortest_new[o, d]

        if (isfinite(sold) != isfinite(snew))
            push!(A, (o, d))
            continue
        end

        if !isfinite(sold) && !isfinite(snew)
            continue
        end

        if abs(sold - snew) > tol
            push!(A, (o, d))
            continue
        end

        stops_old, lines_old = _reconstruct_from_fw(o, d, path_old, lineinfo_old)
        stops_new, lines_new = _reconstruct_from_fw(o, d, path_new, lineinfo_new)

        if stops_old != stops_new
            push!(A, (o, d))
            continue
        end

        if lines_old != lines_new
            push!(A, (o, d))
            continue
        end
    end

    return A
end

# HELPERS
# Add transfer penalty
function add_transfer_penalty(directlink::AbstractMatrix{Float64}, tp::Number)
    n = size(directlink,1)
    D = copy(directlink)
    @inbounds for i in 1:n, j in 1:n
        if i != j && isfinite(D[i,j])
            D[i,j] += tp
        end
    end
    return D
end

# Subtract first boarding penalty
function subtract_first_boarding!(shortest::AbstractMatrix{Float64}, tp::Number)
    n = size(shortest,1)
    @inbounds for i in 1:n, j in 1:n
        if i != j && isfinite(shortest[i,j])
            shortest[i,j] -= tp
        end
    end
    return shortest
end

# FULL EVALUATION 
# Compute full APSP after a move, with split timing.
function full_evaluation_split(network, lineplan_new, tp)

    directlink_new = Matrix{Float64}(undef, 0, 0)
    lineinfo_new   = Matrix{Int}(undef, 0, 0)

    t_dln = @elapsed begin
        directlink_new, lineinfo_new = directlinknetworkinit(network, lineplan_new)
    end

    shortest_new = Matrix{Float64}(undef, size(directlink_new, 1), size(directlink_new, 2))
    path_new     = Matrix{Int}(undef, size(directlink_new, 1), size(directlink_new, 2))

    t_fw_full = @elapsed begin
        dl_tp = add_transfer_penalty(directlink_new, tp)
        shortest_new, path_new = floyd_warshall(dl_tp)
        subtract_first_boarding!(shortest_new, tp)
    end

    return shortest_new, directlink_new, lineinfo_new, path_new, t_dln, t_fw_full
end

# Compute full APSP after a move, returning total time.
function full_evaluation(network, lineplan_new, tp)
    shortest_new, directlink_new, lineinfo_new, path_new, t_dln, t_fw_full =
        full_evaluation_split(network, lineplan_new, tp)
    return shortest_new, directlink_new, lineinfo_new, path_new, (t_dln + t_fw_full)
end

# PARTIAL EVALUATION - DIJKSTRA
# Group OD pairs by origin
function group_by_origin(pairs::Vector{Tuple{Int,Int}})
    G = Dict{Int,Vector{Int}}()
    @inbounds for (o,d) in pairs
        dests = get!(G, o, Int[])
        push!(dests, d)
    end
    return G
end

# Basic Dijkstra
function dijkstra_matrix(cost::AbstractMatrix{Float64}, source::Int)
    n = size(cost,1)
    dist = fill(Inf, n)
    visited = falses(n)

    dist[source] = 0.0

    for _ in 1:n
        u = 0
        best = Inf
        @inbounds for i in 1:n
            if !visited[i] && dist[i] < best
                best = dist[i]
                u = i
            end
        end
        u == 0 && break
        visited[u] = true

        @inbounds for v in 1:n
            if isfinite(cost[u,v])
                alt = dist[u] + cost[u,v]
                if alt < dist[v]
                    dist[v] = alt
                end
            end
        end
    end

    return dist
end

# Local evaluation with basic Dijkstra on detected OD pairs.
function local_dijkstra_evaluation(network,
                                   lineplan_new,
                                   shortest_old,
                                   tp,
                                   A_rule::Set{Tuple{Int,Int}})

    directlink_new = Matrix{Float64}(undef, 0, 0)
    t_dln = @elapsed begin
        directlink_new, _ = directlinknetworkinit(network, lineplan_new)
    end

    pairs = Vector{Tuple{Int,Int}}(collect(A_rule))
    grouped = group_by_origin(pairs)

    dl_tp = add_transfer_penalty(directlink_new, tp)
    shortest_local = copy(shortest_old)

    t_dijkstra = @elapsed begin
        for (o, dests) in grouped
            dist = dijkstra_matrix(dl_tp, o)
            @inbounds for d in dests
                shortest_local[o,d] = dist[d] - tp
            end
        end
    end

    return shortest_local, directlink_new, t_dln, t_dijkstra
end

# Priority-queue Dijkstra using an edge list
function dijkstra_priority_queue(edges_by_origin::Vector{Vector{Tuple{Int,Float64}}}, source::Int)
    n = length(edges_by_origin)
    dist = fill(Inf, n)
    pred = fill(0, n)
    visited = falses(n)

    pq = PriorityQueue{Int,Float64}()

    dist[source] = 0.0
    pq[source] = 0.0

    while !isempty(pq)
        u = dequeue!(pq)

        if visited[u]
            continue
        end
        visited[u] = true

        @inbounds for (v, w) in edges_by_origin[u]
            if !visited[v]
                alt = dist[u] + w
                if alt < dist[v]
                    dist[v] = alt
                    pred[v] = u
                    pq[v] = alt
                end
            end
        end
    end

    return dist, pred
end

# Convert a cost matrix to an edge list grouped by origin
function costmatrix_to_edge_list(cost::AbstractMatrix{Float64})
    n = size(cost, 1)
    edges_by_origin = [Vector{Tuple{Int,Float64}}() for _ in 1:n]

    @inbounds for u in 1:n
        for v in 1:n
            if isfinite(cost[u, v])
                push!(edges_by_origin[u], (v, cost[u, v]))
            end
        end
    end

    return edges_by_origin
end

# Local evaluation with priority-queue Dijkstra on detected OD pairs.
function local_dijkstra_pq_evaluation(network,
                                      lineplan_new,
                                      shortest_old,
                                      tp,
                                      A_rule::Set{Tuple{Int,Int}})

    directlink_new = Matrix{Float64}(undef, 0, 0)
    t_dln = @elapsed begin
        directlink_new, _ = directlinknetworkinit(network, lineplan_new)
    end

    pairs = Vector{Tuple{Int,Int}}(collect(A_rule))
    grouped = group_by_origin(pairs)

    dl_tp = add_transfer_penalty(directlink_new, tp)
    edges_by_origin = costmatrix_to_edge_list(dl_tp)

    shortest_local = copy(shortest_old)

    t_dijkstra_pq = @elapsed begin
        for (o, dests) in grouped
            dist, _ = dijkstra_priority_queue(edges_by_origin, o)
            @inbounds for d in dests
                shortest_local[o, d] = dist[d] - tp
            end
        end
    end

    return shortest_local, directlink_new, t_dln, t_dijkstra_pq
end

# PARTIAL EVALUATION - SUBSET-BASED FLOYD-WARSHALL

const INCLUDE_LINE_STOPS_IN_SUBSET = true

# Extract node subset from OD pairs
function nodes_from_od_pairs(A_rule::Set{Tuple{Int,Int}})
    S = Set{Int}()
    for (o,d) in A_rule
        push!(S,o)
        push!(S,d)
    end
    return sort!(collect(S))
end

# Extract stops for a given subset of lines
function line_stops_for_subset_local(lineplan_old,
                                     lineplan_new,
                                     move_lines)

    nodes = Set{Int}()

    function add_line_if_exists!(lp, line_id)
        if 1 <= line_id <= size(lp, 1)
            stops = Int.(filter(!=(0), lp[line_id, :]))
            for v in stops
                push!(nodes, v)
            end
        end
    end

    if move_lines === nothing
        return nodes

    elseif move_lines isa Int
        add_line_if_exists!(lineplan_old, move_lines)
        add_line_if_exists!(lineplan_new, move_lines)

    elseif move_lines isa Tuple || move_lines isa Vector
        for line_id in move_lines
            add_line_if_exists!(lineplan_old, line_id)
            add_line_if_exists!(lineplan_new, line_id)
        end

    else
        error("Unsupported move_lines type: $(typeof(move_lines))")
    end

    return nodes
end

# Extract node subset from affected OD pairs + line stops
function subset_nodes_from_affected_local(lineplan_old,
                                          lineplan_new,
                                          move_lines,
                                          A_rule::Set{Tuple{Int,Int}},
                                          directlink_new::AbstractMatrix{Float64})

    nodes = Set{Int}()

    for (o, d) in A_rule
        push!(nodes, o)
        push!(nodes, d)
    end

    if INCLUDE_LINE_STOPS_IN_SUBSET
        extra_nodes = line_stops_for_subset_local(lineplan_old, lineplan_new, move_lines)
        union!(nodes, extra_nodes)
    end

    nodes_vec = sort!(collect(nodes))

    return nodes_vec
end

# Local evaluation with subset-based Floyd-Warshall on detected OD pairs.
function local_fw_subset_evaluation(network,
                                    lineplan_old,
                                    lineplan_new,
                                    move_lines,
                                    tp,
                                    A_rule::Set{Tuple{Int,Int}})

    directlink_new = Matrix{Float64}(undef, 0, 0)
    t_dln = @elapsed begin
        directlink_new, _ = directlinknetworkinit(network, lineplan_new)
    end

    nodes = subset_nodes_from_affected_local(
        lineplan_old,
        lineplan_new,
        move_lines,
        A_rule,
        directlink_new
    )

    if length(nodes) < 2
        return nothing, nodes, directlink_new, t_dln, 0.0
    end

    shortest_sub = Matrix{Float64}(undef, length(nodes), length(nodes))
    t_fw_subset = @elapsed begin
        dl_tp = add_transfer_penalty(directlink_new, tp)
        sub_cost = Matrix(dl_tp[nodes, nodes])
        shortest_sub, _ = floyd_warshall(sub_cost)
        subtract_first_boarding!(shortest_sub, tp)
    end

    return shortest_sub, nodes, directlink_new, t_dln, t_fw_subset
end

# RESULTS

# Compare rule vs ground truth
function detection_quality(A_true, A_rule)
    missed = setdiff(A_true, A_rule)
    extra  = setdiff(A_rule, A_true)
    return length(A_true), length(A_rule), length(missed), length(extra)
end

# Check subset FW quality against full FW on predicted pairs
function subset_quality(shortest_full,
                        shortest_sub,
                        nodes_subset,
                        A_rule;
                        tol=1e-9)

    if shortest_sub === nothing
        return 0, 0
    end

    idx = Dict(nodes_subset[i] => i for i in eachindex(nodes_subset))

    bad = 0
    total = 0

    for (o,d) in A_rule
        total += 1
        io = get(idx,o,0)
        id = get(idx,d,0)
        if io==0 || id==0
            continue
        end
        if abs(shortest_sub[io,id] - shortest_full[o,d]) > tol
            bad += 1
        end
    end

    return bad, total
end