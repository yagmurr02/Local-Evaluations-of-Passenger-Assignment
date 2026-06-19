"""
Unified evaluation module supporting both SA and VNS moves.
This file contains all common evaluation algorithms and move detection logic.
"""

include(joinpath(@__DIR__, "common_lineplan_utils.jl"))
include(joinpath(@__DIR__, "..", "evaluator.jl"))
include(joinpath(@__DIR__, "..", "detection", "Rules_StraightforwardApproach.jl"))

using DataStructures: PriorityQueue, dequeue!

# Subset FW configuration
const INCLUDE_LINE_STOPS_IN_SUBSET = true

# OBJECTIVE AND SOLUTION HELPERS

# Return the nodes served by at least one finite shortest-path connection
function get_served_stops(shortest::AbstractMatrix{Float64})
    n_nodes = size(shortest, 1)
    return vec([any(isfinite.(shortest[i, :])) || any(isfinite.(shortest[:, i])) for i in 1:n_nodes])
end

# Create a TransitSolution object from evaluation outputs
function make_solution(lineplan, obj, shortest, directlink, path_old, lineinfo)
    return TransitSolution(
        deepcopy(lineplan),
        obj,
        shortest,
        directlink,
        path_old,
        lineinfo
    )
end

# Compute objective values and shortest-path outputs for a lineplan
function compute_objective_from_shortest(lineplan,
                                         network,
                                         demand,
                                         tp;
                                         min_nodes::Int = 2,
                                         max_nodes::Int = size(lineplan, 2))

    n_nodes = size(network, 1)

    feasible = is_feasible(lineplan, network, demand, n_nodes;
                           tp = tp,
                           min_nodes = min_nodes,
                           max_nodes = max_nodes,
                           enforce_connected = false,
                           enforce_all_covered = true)

    shortest, lineinfo, directlink, path = shortest_paths(network, lineplan, tp)
    served_stops = get_served_stops(shortest)

    if !feasible
        TTT = Inf
        AVT = Inf
    else
        TTT = 0.0

        @inbounds for o in 1:n_nodes, d in 1:n_nodes
            o == d && continue

            dem = demand[o, d]
            dem == 0 && continue

            val = shortest[o, d]

            if !isfinite(val)
                TTT = Inf
                break
            end

            TTT += dem * val
        end

        total_demand = sum(demand)
        AVT = isfinite(TTT) ? TTT / total_demand : Inf
    end

    return TTT, AVT, served_stops, shortest, lineinfo, directlink, path
end

# Compute objective values and shortest-path outputs for a feasible lineplan
function compute_objective_from_shortest_feasible(lineplan,
                                                  network,
                                                  demand,
                                                  tp)

    n_nodes = size(network, 1)

    shortest, lineinfo, directlink, path = shortest_paths(network, lineplan, tp)
    served_stops = get_served_stops(shortest)

    TTT = 0.0

    @inbounds for o in 1:n_nodes, d in 1:n_nodes
        o == d && continue

        dem = demand[o, d]
        dem == 0 && continue

        val = shortest[o, d]

        if !isfinite(val)
            TTT = Inf
            break
        end

        TTT += dem * val
    end

    total_demand = sum(demand)
    AVT = isfinite(TTT) ? TTT / total_demand : Inf

    return TTT, AVT, served_stops, shortest, lineinfo, directlink, path
end

# Perform full evaluation of a lineplan and return objective and shortest-path outputs
function evaluate_full_solution(network, lineplan, demand, tp)
    obj, _, _, shortest, lineinfo, directlink, path =
        compute_objective_from_shortest(lineplan, network, demand, tp;
                                        min_nodes = 2,
                                        max_nodes = size(lineplan, 2))

    return obj, shortest, directlink, path, lineinfo
end

# Perform full evaluation of a feasible lineplan and return objective and shortest-path outputs
function evaluate_full_solution_feasible(network, lineplan, demand, tp)
    obj, _, _, shortest, lineinfo, directlink, path =
        compute_objective_from_shortest_feasible(lineplan, network, demand, tp)

    return obj, shortest, directlink, path, lineinfo
end

# Build a TransitSolution object using full evaluation
function build_solution(lineplan, network, demand, tp)
    obj, shortest, directlink, path_old, lineinfo =
        evaluate_full_solution(network, lineplan, demand, tp)

    return make_solution(lineplan, obj, shortest, directlink, path_old, lineinfo)
end

# Build a TransitSolution object for a feasible lineplan using full evaluation
function build_solution_feasible(lineplan, network, demand, tp)
    obj, shortest, directlink, path_old, lineinfo =
        evaluate_full_solution_feasible(network, lineplan, demand, tp)

    return make_solution(lineplan, obj, shortest, directlink, path_old, lineinfo)
end

# MOVE DETECTION AND PREPROCESSING HELPERS

# Extract the node removed from a terminal-removal move
function removed_terminal_node(lineplan_old, line_id::Int; position::String = "start")
    route = route_from_lineplan(lineplan_old, line_id)
    isempty(route) && error("Line $line_id is empty; cannot determine removed terminal.")

    if position == "start"
        return route[1]
    elseif position == "end"
        return route[end]
    else
        error("Unknown terminal position: $position")
    end
end

# Build transfer-penalized direct-link and line-assignment matrices
function build_directlink_tp(network, lineplan, tp)
    directlink, lineinfo = directlinknetworkinit(network, lineplan)
    n = size(directlink, 1)

    W_tp = copy(directlink)
    @inbounds for i in 1:n, j in 1:n
        if i != j && isfinite(W_tp[i, j])
            W_tp[i, j] += tp
        end
    end

    return directlink, lineinfo, W_tp
end

# Apply the appropriate move-specific detection rule to identify affected OD pairs
function detect_affected_ods(current_sol::TransitSolution,
                             cand_lineplan,
                             move,
                             network,
                             tp)

    lineplan_old   = current_sol.lineplan
    shortest_old   = current_sol.shortest
    directlink_old = current_sol.directlink
    path_old       = current_sol.path_old
    lineinfo_old   = current_sol.lineinfo

    directlink_new, lineinfo_new = directlinknetworkinit(network, cand_lineplan)

    if move isa TerminalInsertion
        return detect_insertion_terminal(
            shortest_old,
            directlink_old,
            directlink_new,
            lineplan_old,
            move.line,
            move.node,
            tp;
            path_old = path_old,
            lineinfo_old = lineinfo_old,
            lineinfo_new = lineinfo_new
        )

    elseif move isa TerminalRemoval
        y = removed_terminal_node(lineplan_old, move.line)

        return detect_removal_terminal(
            shortest_old,
            directlink_old,
            directlink_new,
            lineplan_old,
            move.line,
            y,
            tp;
            path_old = path_old,
            lineinfo_old = lineinfo_old,
            lineinfo_new = lineinfo_new
        )

    elseif move isa Reversal
        return detect_reversal(
            shortest_old,
            directlink_old,
            directlink_new,
            lineplan_old,
            move.line,
            tp;
            path_old = path_old,
            lineinfo_old = lineinfo_old,
            lineinfo_new = lineinfo_new
        )

    elseif move isa VNSInsertion
        return detect_insertion_internal(
            shortest_old,
            directlink_old,
            directlink_new,
            lineplan_old,
            move.line,
            move.left_node,
            move.right_node,
            move.inserted_node,
            tp;
            path_old = path_old,
            lineinfo_old = lineinfo_old,
            lineinfo_new = lineinfo_new
        )

    elseif move isa VNSInsertionTerminal
        return detect_insertion_terminal(
            shortest_old,
            directlink_old,
            directlink_new,
            lineplan_old,
            move.line,
            move.node,
            tp;
            path_old = path_old,
            lineinfo_old = lineinfo_old,
            lineinfo_new = lineinfo_new
        )

    elseif move isa VNSRemoval
        return detect_removal_internal(
            shortest_old,
            directlink_old,
            directlink_new,
            lineplan_old,
            move.line,
            move.removed_node,
            tp;
            path_old = path_old,
            lineinfo_old = lineinfo_old,
            lineinfo_new = lineinfo_new
        )

    elseif move isa VNSRemovalTerminal
        y = removed_terminal_node(lineplan_old, move.line; position = move.position)

        return detect_removal_terminal(
            shortest_old,
            directlink_old,
            directlink_new,
            lineplan_old,
            move.line,
            y,
            tp;
            path_old = path_old,
            lineinfo_old = lineinfo_old,
            lineinfo_new = lineinfo_new
        )

    elseif move isa VNSSwapWithinLine
        return detect_swap_within_line(
            shortest_old,
            directlink_old,
            directlink_new,
            lineplan_old,
            move.line,
            move.node1,
            move.node2,
            tp;
            path_old = path_old,
            lineinfo_old = lineinfo_old,
            lineinfo_new = lineinfo_new
        )

    elseif move isa VNSSubstitution
        return detect_substitution(
            shortest_old,
            directlink_old,
            directlink_new,
            lineplan_old,
            move.line,
            move.old_node,
            move.new_node,
            tp;
            path_old = path_old,
            lineinfo_old = lineinfo_old,
            lineinfo_new = lineinfo_new
        )

    elseif move isa VNSPartialSegmentSwap
        return detect_partial_segment_swap(
            shortest_old,
            directlink_old,
            directlink_new,
            lineplan_old,
            move.line1,
            move.node1a,
            move.node1b,
            move.line2,
            move.node2a,
            move.node2b,
            tp;
            path_old = path_old,
            lineinfo_old = lineinfo_old,
            lineinfo_new = lineinfo_new
        )

    elseif move isa VNSReversal
        return detect_reversal(
            shortest_old,
            directlink_old,
            directlink_new,
            lineplan_old,
            move.line,
            tp;
            path_old = path_old,
            lineinfo_old = lineinfo_old,
            lineinfo_new = lineinfo_new
        )

    elseif move === nothing
        return Set{Tuple{Int,Int}}()

    else
        error("No detection rule defined for move type $(typeof(move))")
    end
end

# Group affected OD pairs by origin node
function affected_pairs_by_origin(affected::Set{Tuple{Int,Int}})
    by_origin = Dict{Int, Vector{Int}}()

    for (o, d) in affected
        if !haskey(by_origin, o)
            by_origin[o] = Int[]
        end
        push!(by_origin[o], d)
    end

    return by_origin
end

# Check whether a periodic full evaluation should be triggered
function should_do_periodic_full_eval(iteration::Int, full_eval_every::Int)
    return full_eval_every > 0 && iteration > 0 && (iteration % full_eval_every == 0)
end

# PATH RECONSTRUCTION HELPERS

# Reconstruct a path from a Dijkstra parent vector
function reconstruct_path_from_parent(parent::Vector{Int}, src::Int, dst::Int)
    src == dst && return [src]

    if parent[dst] == 0
        return Int[]
    end

    path = Int[]
    v = dst

    while v != 0 && v != src
        push!(path, v)
        v = parent[v]
    end

    if v == 0
        return Int[]
    end

    push!(path, src)
    reverse!(path)

    return path
end

# Encode a stop sequence in Floyd–Warshall intermediate-node format
function encode_path_as_fw_split!(pathmat::AbstractMatrix{Int}, stops::Vector{Int})
    n = length(stops)

    if n <= 1
        return
    elseif n == 2
        pathmat[stops[1], stops[2]] = 0
        return
    else
        mid = fld(n + 1, 2)
        i = stops[1]
        j = stops[end]
        k = stops[mid]

        pathmat[i, j] = k

        encode_path_as_fw_split!(pathmat, stops[1:mid])
        encode_path_as_fw_split!(pathmat, stops[mid:end])

        return
    end
end

# Update the Floyd–Warshall path matrix for one OD pair
function update_fw_path_for_od!(pathmat::AbstractMatrix{Int},
                                stops::Vector{Int},
                                o::Int,
                                d::Int)
    if isempty(stops)
        pathmat[o, d] = 0
        return
    elseif length(stops) == 1
        pathmat[o, d] = 0
        return
    else
        encode_path_as_fw_split!(pathmat, stops)
        return
    end
end

# Reconstruct a path from subset Floyd–Warshall predecessor information
function reconstruct_path_from_subset_pred(pred_sub::AbstractMatrix{Int},
                                           nodes::Vector{Int},
                                           io::Int,
                                           id::Int)
    io == id && return [nodes[io]]

    if pred_sub[io, id] == 0
        return Int[]
    end

    path_idx = Int[]
    v = id

    while v != 0 && v != io
        push!(path_idx, v)
        v = pred_sub[io, v]
    end

    if v == 0
        return Int[]
    end

    push!(path_idx, io)
    reverse!(path_idx)

    return [nodes[t] for t in path_idx]
end

# Reconstruct the full stop sequence from a Floyd–Warshall path matrix
function reconstruct_fw_path_nodes(pathmat::AbstractMatrix{Int},
                                   shortest::AbstractMatrix{Float64},
                                   o::Int,
                                   d::Int)

    if o == d
        return [o]
    end

    if !isfinite(shortest[o, d])
        return Int[]
    end

    k = pathmat[o, d]

    if k == 0
        return [o, d]
    else
        left = reconstruct_fw_path_nodes(pathmat, shortest, o, k)
        right = reconstruct_fw_path_nodes(pathmat, shortest, k, d)

        isempty(left)  && return Int[]
        isempty(right) && return Int[]

        return vcat(left[1:end-1], right)
    end
end

# DIJKSTRA HELPERS

# Build an adjacency-list representation from a cost matrix
function build_edges_by_origin(cost::AbstractMatrix{Float64})
    n = size(cost, 1)
    edges_by_origin = [Vector{Tuple{Int,Float64}}() for _ in 1:n]

    @inbounds for u in 1:n
        for v in 1:n
            w = cost[u, v]
            if u != v && isfinite(w)
                push!(edges_by_origin[u], (v, w))
            end
        end
    end

    return edges_by_origin
end

# Compute shortest paths using matrix-based Dijkstra with transfer-penalized costs
function dijkstra_basic_tp(cost::AbstractMatrix{Float64}, source::Int; tol::Float64 = 0.0)
    n = size(cost, 1)
    dist = fill(Inf, n)
    parent = zeros(Int, n)
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
            if isfinite(cost[u, v])
                alt = dist[u] + cost[u, v]
                if alt + tol < dist[v]
                    dist[v] = alt
                    parent[v] = u
                end
            end
        end
    end

    return dist, parent
end

# Compute shortest paths using priority-queue Dijkstra with transfer-penalized costs
function dijkstra_pq_tp(edges_by_origin::Vector{Vector{Tuple{Int,Float64}}}, source::Int; tol::Float64 = 0.0)
    n = length(edges_by_origin)
    dist = fill(Inf, n)
    parent = fill(0, n)
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
                if alt + tol < dist[v]
                    dist[v] = alt
                    parent[v] = u
                    pq[v] = alt
                end
            end
        end
    end

    return dist, parent
end

# SUBSET FLOYD–WARSHALL HELPERS

# Collect line stops affected by the move for subset construction
function line_stops_for_subset(lineplan_old, cand_lineplan, move)
    nodes = Set{Int}()

    if move === nothing
        return nodes
    end

    if (move isa TerminalInsertion || move isa TerminalRemoval || move isa Reversal ||
        move isa VNSInsertion ||
        move isa VNSInsertionTerminal ||
        move isa VNSRemoval ||
        move isa VNSRemovalTerminal ||
        move isa VNSSwapWithinLine ||
        move isa VNSSubstitution ||
        move isa VNSReversal)

        line_id = move.line

        old_route = route_from_lineplan(lineplan_old, line_id)
        new_route = route_from_lineplan(cand_lineplan, line_id)

        for v in old_route
            push!(nodes, v)
        end
        for v in new_route
            push!(nodes, v)
        end

    elseif move isa VNSPartialSegmentSwap
        for line_id in (move.line1, move.line2)
            old_route = route_from_lineplan(lineplan_old, line_id)
            new_route = route_from_lineplan(cand_lineplan, line_id)

            for v in old_route
                push!(nodes, v)
            end
            for v in new_route
                push!(nodes, v)
            end
        end
    end

    return nodes
end

# Build the node subset from affected OD endpoints and optionally affected line stops
function subset_nodes_from_affected(current_sol::TransitSolution,
                                    cand_lineplan,
                                    move,
                                    affected::Set{Tuple{Int,Int}},
                                    network,
                                    tp;
                                    expand_one_hop::Bool = false)
    nodes = Set{Int}()

    for (o, d) in affected
        push!(nodes, o)
        push!(nodes, d)
    end

    if INCLUDE_LINE_STOPS_IN_SUBSET
        extra_nodes = line_stops_for_subset(current_sol.lineplan, cand_lineplan, move)
        union!(nodes, extra_nodes)
    end

    if expand_one_hop && !isempty(nodes)
        directlink_new, _, _ = build_directlink_tp(network, cand_lineplan, tp)
        n = size(directlink_new, 1)

        one_hop_seeds = collect(nodes)
        one_hop_nodes = Set{Int}()

        for v in one_hop_seeds
            @inbounds for u in 1:n
                if u != v && (isfinite(directlink_new[v, u]) || isfinite(directlink_new[u, v]))
                    push!(one_hop_nodes, u)
                end
            end
        end

        union!(nodes, one_hop_nodes)
    end

    return sort!(collect(nodes))
end

# Extract the induced cost submatrix for the selected node subset
function induced_submatrix(W::AbstractMatrix{Float64}, nodes::Vector{Int})
    k = length(nodes)
    Wsub = fill(Inf, k, k)

    @inbounds for i in 1:k
        Wsub[i, i] = 0.0
        gi = nodes[i]
        for j in 1:k
            gj = nodes[j]
            Wsub[i, j] = W[gi, gj]
        end
    end

    return Wsub
end

# Count the number of unique origins among affected OD pairs
function count_affected_origins(affected::Set{Tuple{Int,Int}})
    origins = Set{Int}()
    for (o, _) in affected
        push!(origins, o)
    end
    return length(origins)
end

# Count finite directed arcs in a cost matrix
function count_finite_arcs(W::AbstractMatrix{Float64})
    n = size(W, 1)
    m = 0

    @inbounds for i in 1:n, j in 1:n
        if i != j && isfinite(W[i, j])
            m += 1
        end
    end

    return m
end

# Estimate whether Dijkstra is cheaper than subset-based Floyd–Warshall
function dijkstra_cheaper_than_subset_fw(affected::Set{Tuple{Int,Int}},
                                         current_sol::TransitSolution,
                                         cand_lineplan,
                                         move,
                                         network,
                                         tp;
                                         dijkstra_mode::Symbol = :pq)

    isempty(affected) && return true

    n = size(network, 1)
    r = count_affected_origins(affected)
    subset_nodes = subset_nodes_from_affected(current_sol, cand_lineplan, move, affected, network, tp)
    s = length(subset_nodes)

    s == 0 && return true

    subset_fw_work = float(s)^3

    if dijkstra_mode == :basic
        dijkstra_work = float(r) * float(n)^2

    elseif dijkstra_mode == :pq
        _, _, W_tp = build_directlink_tp(network, cand_lineplan, tp)
        m = count_finite_arcs(W_tp)
        logn = max(log(float(max(n, 2))), 1.0)
        dijkstra_work = float(r) * float(m) * logn

    else
        error("Unknown dijkstra_mode: $dijkstra_mode")
    end

    return dijkstra_work <= subset_fw_work
end

# Compute shortest paths on a node subset using Floyd–Warshall with predecessor tracking
function floyd_warshall_subset_tp(Wsub::AbstractMatrix{Float64}; tol::Float64 = 0.0)
    k = size(Wsub, 1)

    dist = copy(Wsub)
    pred = zeros(Int, k, k)

    @inbounds for i in 1:k, j in 1:k
        if i != j && isfinite(dist[i, j])
            pred[i, j] = i
        end
    end

    @inbounds for m in 1:k
        for i in 1:k
            dim = dist[i, m]
            isfinite(dim) || continue

            for j in 1:k
                dmj = dist[m, j]
                isfinite(dmj) || continue

                nd = dim + dmj
                if nd + tol < dist[i, j]
                    dist[i, j] = nd
                    pred[i, j] = pred[m, j]
                end
            end
        end
    end

    return dist, pred
end

# CHANGED-ARC AND PATH ANALYSIS HELPERS

# Identify direct links whose existence, travel time, or line assignment changed
function changed_direct_arcs(directlink_old::AbstractMatrix{Float64},
                             directlink_new::AbstractMatrix{Float64},
                             lineinfo_old,
                             lineinfo_new;
                             tol::Float64 = 1e-9)

    n = size(directlink_old, 1)
    changed = falses(n, n)

    @inbounds for i in 1:n, j in 1:n
        oldv = directlink_old[i, j]
        newv = directlink_new[i, j]

        # Direct-link existence changed
        if isfinite(oldv) != isfinite(newv)
            changed[i, j] = true
            continue
        end

        # Direct-link time changed
        if isfinite(oldv) && isfinite(newv) && abs(oldv - newv) > tol
            changed[i, j] = true
            continue
        end

        # Same direct-link time/existence, but line assignment changed
        if isfinite(oldv) && isfinite(newv)
            if lineinfo_old[i, j] != lineinfo_new[i, j]
                changed[i, j] = true
            end
        end
    end

    return changed
end

# Check whether a path uses at least one changed direct arc
function path_uses_changed_arc(path_nodes::Vector{Int},
                               changed_arcs::AbstractMatrix{Bool})
    length(path_nodes) <= 1 && return false

    @inbounds for t in 1:(length(path_nodes) - 1)
        u = path_nodes[t]
        v = path_nodes[t + 1]
        if changed_arcs[u, v]
            return true
        end
    end

    return false
end

# SHARED OD UPDATE HELPERS

# Update one OD pair in the local solution, including shortest value, path encoding, and objective contribution
function apply_od_update!(shortest_new,
                          path_new,
                          obj_new_ref,
                          current_sol::TransitSolution,
                          demand,
                          o::Int,
                          d::Int,
                          new_val,
                          path_nodes::Vector{Int})

    old_val = current_sol.shortest[o, d]
    shortest_new[o, d] = new_val

    if !isfinite(new_val) || o == d
        path_new[o, d] = 0
    else
        update_fw_path_for_od!(path_new, path_nodes, o, d)
    end

    if o != d
        dem = demand[o, d]

        if dem > 0
            if !isfinite(new_val)
                obj_new_ref[] = Inf
                return
            end

            if isfinite(obj_new_ref[])
                if isfinite(old_val)
                    obj_new_ref[] -= dem * old_val
                end
                obj_new_ref[] += dem * new_val
            end
        end
    end
end

# LOCAL EVALUATION METHODS

# Evaluate affected OD pairs using local Dijkstra and update the solution outputs
function evaluate_local_dijkstra_from_affected(current_sol::TransitSolution,
                                               cand_lineplan,
                                               affected::Set{Tuple{Int,Int}},
                                               network,
                                               demand,
                                               tp;
                                               dijkstra_mode::Symbol = :pq)

    directlink_new, lineinfo_new, W_tp = build_directlink_tp(network, cand_lineplan, tp)

    edges_by_origin =
        dijkstra_mode == :pq ? build_edges_by_origin(W_tp) : Vector{Vector{Tuple{Int,Float64}}}()

    shortest_new = copy(current_sol.shortest)
    path_new     = copy(current_sol.path_old)
    obj_new_ref  = Ref(current_sol.obj)

    if !isempty(affected)
        by_origin = affected_pairs_by_origin(affected)

        for (o, dests) in by_origin
            if dijkstra_mode == :basic
                dist_tp, parent = dijkstra_basic_tp(W_tp, o)
            elseif dijkstra_mode == :pq
                dist_tp, parent = dijkstra_pq_tp(edges_by_origin, o)
            else
                error("Unknown dijkstra_mode: $dijkstra_mode")
            end

            for d in dests
                new_tp  = dist_tp[d]
                new_val = (o == d || !isfinite(new_tp)) ? new_tp : (new_tp - tp)
                stops = (!isfinite(new_val) || o == d) ? Int[] : reconstruct_path_from_parent(parent, o, d)

                apply_od_update!(shortest_new, path_new, obj_new_ref, current_sol, demand, o, d, new_val, stops)
            end
        end
    end

    return obj_new_ref[], shortest_new, directlink_new, path_new, lineinfo_new
end

# Evaluate affected OD pairs using subset-based Floyd–Warshall and update the solution outputs
function evaluate_local_floyd_subset_from_affected(current_sol::TransitSolution,
                                                   cand_lineplan,
                                                   move,
                                                   affected::Set{Tuple{Int,Int}},
                                                   network,
                                                   demand,
                                                   tp)

    directlink_new, lineinfo_new, W_tp = build_directlink_tp(network, cand_lineplan, tp)

    shortest_new = copy(current_sol.shortest)
    path_new     = copy(current_sol.path_old)
    obj_new_ref  = Ref(current_sol.obj)

    isempty(affected) && return obj_new_ref[], shortest_new, directlink_new, path_new, lineinfo_new

    nodes = subset_nodes_from_affected(current_sol, cand_lineplan, move, affected, network, tp)
    node_to_sub = Dict(node => idx for (idx, node) in enumerate(nodes))

    Wsub = induced_submatrix(W_tp, nodes)
    dist_sub, pred_sub = floyd_warshall_subset_tp(Wsub)

    for (o, d) in affected
        io = node_to_sub[o]
        id = node_to_sub[d]

        new_tp  = dist_sub[io, id]
        new_val = (o == d || !isfinite(new_tp)) ? new_tp : (new_tp - tp)
        stops = (!isfinite(new_val) || o == d) ? Int[] : reconstruct_path_from_subset_pred(pred_sub, nodes, io, id)

        apply_od_update!(shortest_new, path_new, obj_new_ref, current_sol, demand, o, d, new_val, stops)
    end

    return obj_new_ref[], shortest_new, directlink_new, path_new, lineinfo_new
end

# LOCAL EVALUATION WRAPPERS

# Evaluate a candidate solution using local matrix-based Dijkstra
function evaluate_local_dijkstra_basic_solution(current_sol::TransitSolution,
                                                cand_lineplan,
                                                move,
                                                network,
                                                demand,
                                                tp)

    local affected
    detection_time = @elapsed begin
        affected = detect_affected_ods(current_sol, cand_lineplan, move, network, tp)
    end

    affected === nothing && (affected = Set{Tuple{Int,Int}}())

    obj, shortest, directlink, path_new, lineinfo =
        evaluate_local_dijkstra_from_affected(
            current_sol, cand_lineplan, affected, network, demand, tp;
            dijkstra_mode = :basic
        )

    return obj, shortest, directlink, path_new, lineinfo, detection_time
end

# Build a TransitSolution using local matrix-based Dijkstra evaluation
function build_solution_local_dijkstra_basic(lineplan,
                                             current_sol::TransitSolution,
                                             move,
                                             network,
                                             demand,
                                             tp)

    obj, shortest, directlink, path_old, lineinfo, detection_time =
        evaluate_local_dijkstra_basic_solution(current_sol, lineplan, move, network, demand, tp)

    return make_solution(lineplan, obj, shortest, directlink, path_old, lineinfo), detection_time
end

# Evaluate a candidate solution using local priority-queue Dijkstra
function evaluate_local_dijkstra_pq_solution(current_sol::TransitSolution,
                                             cand_lineplan,
                                             move,
                                             network,
                                             demand,
                                             tp)

    local affected
    detection_time = @elapsed begin
        affected = detect_affected_ods(current_sol, cand_lineplan, move, network, tp)
    end

    affected === nothing && (affected = Set{Tuple{Int,Int}}())

    obj, shortest, directlink, path_new, lineinfo =
        evaluate_local_dijkstra_from_affected(
            current_sol, cand_lineplan, affected, network, demand, tp;
            dijkstra_mode = :pq
        )

    return obj, shortest, directlink, path_new, lineinfo, detection_time
end

# Build a TransitSolution using local priority-queue Dijkstra evaluation
function build_solution_local_dijkstra_pq(lineplan,
                                          current_sol::TransitSolution,
                                          move,
                                          network,
                                          demand,
                                          tp)

    obj, shortest, directlink, path_old, lineinfo, detection_time =
        evaluate_local_dijkstra_pq_solution(current_sol, lineplan, move, network, demand, tp)

    return make_solution(lineplan, obj, shortest, directlink, path_old, lineinfo), detection_time
end

# Evaluate a candidate solution using local subset-based Floyd–Warshall
function evaluate_local_floyd_subset_solution(current_sol::TransitSolution,
                                              cand_lineplan,
                                              move,
                                              network,
                                              demand,
                                              tp)

    local affected
    detection_time = @elapsed begin
        affected = detect_affected_ods(current_sol, cand_lineplan, move, network, tp)
    end

    affected === nothing && (affected = Set{Tuple{Int,Int}}())

    obj, shortest, directlink, path_new, lineinfo =
        evaluate_local_floyd_subset_from_affected(
            current_sol, cand_lineplan, move, affected, network, demand, tp
        )

    return obj, shortest, directlink, path_new, lineinfo, detection_time
end

# Build a TransitSolution using local subset-based Floyd–Warshall evaluation
function build_solution_local_floyd_subset(lineplan,
                                           current_sol::TransitSolution,
                                           move,
                                           network,
                                           demand,
                                           tp)

    obj, shortest, directlink, path_old, lineinfo, detection_time =
        evaluate_local_floyd_subset_solution(current_sol, lineplan, move, network, demand, tp)

    return make_solution(lineplan, obj, shortest, directlink, path_old, lineinfo), detection_time
end

# SUBSET FW PREPROCESSING

# Precompute subset-based Floyd–Warshall data for local evaluation
function compute_subset_fw_data(current_sol::TransitSolution,
                                cand_lineplan,
                                move,
                                affected::Set{Tuple{Int,Int}},
                                network,
                                tp)

    directlink_new, lineinfo_new, W_tp = build_directlink_tp(network, cand_lineplan, tp)

    if isempty(affected)
        return directlink_new, lineinfo_new, Int[], Dict{Int,Int}(), zeros(Float64, 0, 0), zeros(Int, 0, 0)
    end

    nodes = subset_nodes_from_affected(current_sol, cand_lineplan, move, affected, network, tp)
    node_to_sub = Dict(node => idx for (idx, node) in enumerate(nodes))

    Wsub = induced_submatrix(W_tp, nodes)
    dist_sub, pred_sub = floyd_warshall_subset_tp(Wsub)

    return directlink_new, lineinfo_new, nodes, node_to_sub, dist_sub, pred_sub
end

# PATH-INFORMED SUBSET FLOYD–WARSHALL METHODS

# Evaluate a candidate solution using path-informed subset-based Floyd–Warshall
function evaluate_local_floyd_subset_path_solution(current_sol::TransitSolution,
                                                   cand_lineplan,
                                                   move,
                                                   network,
                                                   demand,
                                                   tp)

    local affected
    detection_time = @elapsed begin
        affected = detect_affected_ods(current_sol, cand_lineplan, move, network, tp)
    end

    affected === nothing && (affected = Set{Tuple{Int,Int}}())

    directlink_new, lineinfo_new, nodes, node_to_sub, dist_sub, pred_sub =
        compute_subset_fw_data(current_sol, cand_lineplan, move, affected, network, tp)

    shortest_new = copy(current_sol.shortest)
    path_new     = copy(current_sol.path_old)
    obj_new      = current_sol.obj

    isempty(affected) && return obj_new, shortest_new, directlink_new, path_new, lineinfo_new, detection_time

    changed_arcs = changed_direct_arcs(
        current_sol.directlink,
        directlink_new,
        current_sol.lineinfo,
        lineinfo_new
    )

    for (o, d) in affected
        old_val = current_sol.shortest[o, d]

        io = node_to_sub[o]
        id = node_to_sub[d]

        subset_tp  = dist_sub[io, id]
        subset_val = (o == d || !isfinite(subset_tp)) ? subset_tp : (subset_tp - tp)

        old_path_nodes = reconstruct_fw_path_nodes(current_sol.path_old, current_sol.shortest, o, d)
        old_path_changed = path_uses_changed_arc(old_path_nodes, changed_arcs)

        if !old_path_changed
            if isfinite(old_val) && isfinite(subset_val)
                new_val = min(old_val, subset_val)
            elseif isfinite(old_val)
                new_val = old_val
            else
                new_val = subset_val
            end
        else
            new_val = subset_val
        end

        shortest_new[o, d] = new_val

        if !isfinite(new_val)
            path_new[o, d] = 0
        elseif o == d
            path_new[o, d] = 0
        elseif !old_path_changed && isfinite(old_val) &&
               ((!isfinite(subset_val) && isfinite(old_val)) ||
                (isfinite(old_val) && isfinite(subset_val) && old_val <= subset_val))
            update_fw_path_for_od!(path_new, old_path_nodes, o, d)
        else
            subset_nodes_path = reconstruct_path_from_subset_pred(pred_sub, nodes, io, id)
            update_fw_path_for_od!(path_new, subset_nodes_path, o, d)
        end

        if o != d
            dem = demand[o, d]

            if dem > 0
                if !isfinite(new_val)
                    obj_new = Inf
                elseif isfinite(obj_new)
                    if isfinite(old_val)
                        obj_new -= dem * old_val
                    end
                    obj_new += dem * new_val
                end
            end
        end
    end

    return obj_new, shortest_new, directlink_new, path_new, lineinfo_new, detection_time
end

# Build a TransitSolution using path-informed subset-based Floyd–Warshall evaluation
function build_solution_local_floyd_subset_path(lineplan,
                                                current_sol::TransitSolution,
                                                move,
                                                network,
                                                demand,
                                                tp)

    obj, shortest, directlink, path_old, lineinfo, detection_time =
        evaluate_local_floyd_subset_path_solution(current_sol, lineplan, move, network, demand, tp)

    return make_solution(lineplan, obj, shortest, directlink, path_old, lineinfo), detection_time
end

# Evaluate a candidate solution using path-informed subset-based Floyd–Warshall with Dijkstra fallback
function evaluate_local_floyd_subset_path_fallback_solution(current_sol::TransitSolution,
                                                            cand_lineplan,
                                                            move,
                                                            network,
                                                            demand,
                                                            tp;
                                                            dijkstra_mode::Symbol = :pq)

    local affected
    detection_time = @elapsed begin
        affected = detect_affected_ods(current_sol, cand_lineplan, move, network, tp)
    end

    affected === nothing && (affected = Set{Tuple{Int,Int}}())

    directlink_new, lineinfo_new, nodes, node_to_sub, dist_sub, pred_sub =
        compute_subset_fw_data(current_sol, cand_lineplan, move, affected, network, tp)

    shortest_new = copy(current_sol.shortest)
    path_new     = copy(current_sol.path_old)
    obj_new      = current_sol.obj

    isempty(affected) && return obj_new, shortest_new, directlink_new, path_new, lineinfo_new, detection_time

    changed_arcs = changed_direct_arcs(
        current_sol.directlink,
        directlink_new,
        current_sol.lineinfo,
        lineinfo_new
    )
    
    exact_by_origin = Dict{Int, Vector{Int}}()

    _, _, W_tp = build_directlink_tp(network, cand_lineplan, tp)
    edges_by_origin =
        dijkstra_mode == :pq ? build_edges_by_origin(W_tp) : Vector{Vector{Tuple{Int,Float64}}}()

    for (o, d) in affected
        old_val = current_sol.shortest[o, d]

        io = node_to_sub[o]
        id = node_to_sub[d]

        subset_tp  = dist_sub[io, id]
        subset_val = (o == d || !isfinite(subset_tp)) ? subset_tp : (subset_tp - tp)

        old_path_nodes = reconstruct_fw_path_nodes(current_sol.path_old, current_sol.shortest, o, d)
        old_path_changed = path_uses_changed_arc(old_path_nodes, changed_arcs)

        if !old_path_changed
            if isfinite(old_val) && isfinite(subset_val)
                new_val = min(old_val, subset_val)
            elseif isfinite(old_val)
                new_val = old_val
            else
                new_val = subset_val
            end

            shortest_new[o, d] = new_val

            if !isfinite(new_val)
                path_new[o, d] = 0
            elseif o == d
                path_new[o, d] = 0
            elseif isfinite(old_val) &&
                   ((!isfinite(subset_val) && isfinite(old_val)) ||
                    (isfinite(old_val) && isfinite(subset_val) && old_val <= subset_val))
                update_fw_path_for_od!(path_new, old_path_nodes, o, d)
            else
                subset_nodes_path = reconstruct_path_from_subset_pred(pred_sub, nodes, io, id)
                update_fw_path_for_od!(path_new, subset_nodes_path, o, d)
            end

            if o != d
                dem = demand[o, d]

                if dem > 0
                    if !isfinite(new_val)
                        obj_new = Inf
                    elseif isfinite(obj_new)
                        if isfinite(old_val)
                            obj_new -= dem * old_val
                        end
                        obj_new += dem * new_val
                    end
                end
            end
        else
            if !haskey(exact_by_origin, o)
                exact_by_origin[o] = Int[]
            end
            push!(exact_by_origin[o], d)
        end
    end

    for (o, dests) in exact_by_origin
        if dijkstra_mode == :basic
            dist_tp, parent = dijkstra_basic_tp(W_tp, o)
        elseif dijkstra_mode == :pq
            dist_tp, parent = dijkstra_pq_tp(edges_by_origin, o)
        else
            error("Unknown dijkstra_mode: $dijkstra_mode")
        end

        for d in dests
            old_val = current_sol.shortest[o, d]
            new_tp  = dist_tp[d]
            new_val = (o == d || !isfinite(new_tp)) ? new_tp : (new_tp - tp)

            shortest_new[o, d] = new_val

            if !isfinite(new_val)
                path_new[o, d] = 0
            elseif o == d
                path_new[o, d] = 0
            else
                stops = reconstruct_path_from_parent(parent, o, d)
                update_fw_path_for_od!(path_new, stops, o, d)
            end

            if o != d
                dem = demand[o, d]

                if dem > 0
                    if !isfinite(new_val)
                        obj_new = Inf
                    elseif isfinite(obj_new)
                        if isfinite(old_val)
                            obj_new -= dem * old_val
                        end
                        obj_new += dem * new_val
                    end
                end
            end
        end
    end

    return obj_new, shortest_new, directlink_new, path_new, lineinfo_new, detection_time
end

# Build a TransitSolution using path-informed subset-based Floyd–Warshall with Dijkstra fallback
function build_solution_local_floyd_subset_path_fallback(lineplan,
                                                         current_sol::TransitSolution,
                                                         move,
                                                         network,
                                                         demand,
                                                         tp;
                                                         dijkstra_mode::Symbol = :pq)

    obj, shortest, directlink, path_old, lineinfo, detection_time =
        evaluate_local_floyd_subset_path_fallback_solution(
            current_sol, lineplan, move, network, demand, tp;
            dijkstra_mode = dijkstra_mode
        )

    return make_solution(lineplan, obj, shortest, directlink, path_old, lineinfo), detection_time
end

# Build a TransitSolution using path-informed subset-based Floyd–Warshall with matrix-based Dijkstra fallback
function build_solution_local_fw_subset_path_basic(lineplan,
                                                   current_sol::TransitSolution,
                                                   move,
                                                   network,
                                                   demand,
                                                   tp)

    return build_solution_local_floyd_subset_path_fallback(
        lineplan, current_sol, move, network, demand, tp;
        dijkstra_mode = :basic
    )
end

# Build a TransitSolution using path-informed subset-based Floyd–Warshall with priority-queue Dijkstra fallback
function build_solution_local_fw_subset_path_pq(lineplan,
                                                current_sol::TransitSolution,
                                                move,
                                                network,
                                                demand,
                                                tp)

    return build_solution_local_floyd_subset_path_fallback(
        lineplan, current_sol, move, network, demand, tp;
        dijkstra_mode = :pq
    )
end

# HYBRID LOCAL EVALUATION METHODS

# Evaluate a candidate solution using hybrid local evaluation with Dijkstra fallback
function evaluate_local_hybrid_solution(current_sol::TransitSolution,
                                        cand_lineplan,
                                        move,
                                        network,
                                        demand,
                                        tp;
                                        dijkstra_mode::Symbol = :pq)

    obj_fw, shortest_fw, directlink_fw, path_fw, lineinfo_fw, detection_time_fw =
        evaluate_local_floyd_subset_solution(
            current_sol, cand_lineplan, move, network, demand, tp
        )

    if obj_fw < current_sol.obj
        return obj_fw, shortest_fw, directlink_fw, path_fw, lineinfo_fw, false, detection_time_fw
    end

    if dijkstra_mode == :basic
        obj_dj, shortest_dj, directlink_dj, path_dj, lineinfo_dj, detection_time_dj =
            evaluate_local_dijkstra_basic_solution(
                current_sol, cand_lineplan, move, network, demand, tp
            )
    elseif dijkstra_mode == :pq
        obj_dj, shortest_dj, directlink_dj, path_dj, lineinfo_dj, detection_time_dj =
            evaluate_local_dijkstra_pq_solution(
                current_sol, cand_lineplan, move, network, demand, tp
            )
    else
        error("Unknown dijkstra_mode: $dijkstra_mode")
    end

    return obj_dj, shortest_dj, directlink_dj, path_dj, lineinfo_dj, true, detection_time_fw + detection_time_dj
end

# Build a TransitSolution using hybrid local evaluation with matrix-based Dijkstra fallback
function build_solution_local_hybrid_basic(lineplan,
                                           current_sol::TransitSolution,
                                           move,
                                           network,
                                           demand,
                                           tp)

    obj, shortest, directlink, path_old, lineinfo, used_fallback, detection_time =
        evaluate_local_hybrid_solution(
            current_sol, lineplan, move, network, demand, tp;
            dijkstra_mode = :basic
        )

    return make_solution(lineplan, obj, shortest, directlink, path_old, lineinfo), used_fallback, detection_time
end

# Build a TransitSolution using hybrid local evaluation with priority-queue Dijkstra fallback
function build_solution_local_hybrid_pq(lineplan,
                                        current_sol::TransitSolution,
                                        move,
                                        network,
                                        demand,
                                        tp)

    obj, shortest, directlink, path_old, lineinfo, used_fallback, detection_time =
        evaluate_local_hybrid_solution(
            current_sol, lineplan, move, network, demand, tp;
            dijkstra_mode = :pq
        )

    return make_solution(lineplan, obj, shortest, directlink, path_old, lineinfo), used_fallback, detection_time
end

# THRESHOLD-ADAPTIVE LOCAL EVALUATION

# Evaluate a candidate solution using threshold-adaptive selection among Dijkstra, subset FW, and full evaluation
function evaluate_local_threshold_adaptive_solution(current_sol::TransitSolution,
                                                    cand_lineplan,
                                                    move,
                                                    network,
                                                    demand,
                                                    tp;
                                                    full_threshold::Float64 = 0.70,
                                                    dijkstra_mode::Symbol = :pq)

    if !(0.0 <= full_threshold <= 1.0)
        error("full_threshold must satisfy 0.0 <= full_threshold <= 1.0")
    end

    local affected
    detection_time = @elapsed begin
        affected = detect_affected_ods(current_sol, cand_lineplan, move, network, tp)
    end

    affected === nothing && (affected = Set{Tuple{Int,Int}}())

    total_nodes = size(network, 1)
    subset_nodes = subset_nodes_from_affected(current_sol, cand_lineplan, move, affected, network, tp)
    subset_ratio = total_nodes == 0 ? 0.0 : length(subset_nodes) / total_nodes

    local obj, shortest, directlink, path, lineinfo
    local method_used::Symbol

    if isempty(affected)
        obj        = current_sol.obj
        shortest   = copy(current_sol.shortest)
        directlink = current_sol.directlink
        path       = copy(current_sol.path_old)
        lineinfo   = current_sol.lineinfo
        method_used = :none

    elseif dijkstra_cheaper_than_subset_fw(
               affected, current_sol, cand_lineplan, move, network, tp;
               dijkstra_mode = dijkstra_mode
           )
        obj, shortest, directlink, path, lineinfo =
            evaluate_local_dijkstra_from_affected(
                current_sol, cand_lineplan, affected, network, demand, tp;
                dijkstra_mode = dijkstra_mode
            )
        method_used = :local_dijkstra

    elseif subset_ratio < full_threshold
        obj, shortest, directlink, path, lineinfo =
            evaluate_local_floyd_subset_from_affected(
                current_sol, cand_lineplan, move, affected, network, demand, tp
            )
        method_used = :local_floyd_subset

    else
        obj, shortest, directlink, path, lineinfo =
            evaluate_full_solution_feasible(network, cand_lineplan, demand, tp)
        method_used = :full
    end

    return obj, shortest, directlink, path, lineinfo, detection_time, method_used
end

# Build a TransitSolution using threshold-adaptive local evaluation
function build_solution_local_threshold_adaptive(lineplan,
                                                 current_sol::TransitSolution,
                                                 move,
                                                 network,
                                                 demand,
                                                 tp;
                                                 full_threshold::Float64 = 0.70,
                                                 dijkstra_mode::Symbol = :pq)

    obj, shortest, directlink, path_old, lineinfo, detection_time, method_used =
        evaluate_local_threshold_adaptive_solution(
            current_sol,
            lineplan,
            move,
            network,
            demand,
            tp;
            full_threshold = full_threshold,
            dijkstra_mode = dijkstra_mode
        )

    return make_solution(lineplan, obj, shortest, directlink, path_old, lineinfo), detection_time, method_used
end

# Build a TransitSolution using threshold-adaptive evaluation with matrix-based Dijkstra
function build_solution_local_threshold_basic(lineplan,
                                              current_sol::TransitSolution,
                                              move,
                                              network,
                                              demand,
                                              tp;
                                              full_threshold::Float64 = 0.70)

    return build_solution_local_threshold_adaptive(
        lineplan, current_sol, move, network, demand, tp;
        full_threshold = full_threshold,
        dijkstra_mode = :basic
    )
end

# Build a TransitSolution using threshold-adaptive evaluation with priority-queue Dijkstra
function build_solution_local_threshold_pq(lineplan,
                                           current_sol::TransitSolution,
                                           move,
                                           network,
                                           demand,
                                           tp;
                                           full_threshold::Float64 = 0.70)

    return build_solution_local_threshold_adaptive(
        lineplan, current_sol, move, network, demand, tp;
        full_threshold = full_threshold,
        dijkstra_mode = :pq
    )
end

# For adaptive method testing
# function is_loop_line(lineplan, line::Int)
#     route = route_from_lineplan(lineplan, line)
#     return length(route) != length(unique(route))
# end

# function build_solution_vns_adaptive(lineplan,
#                                      current_sol::TransitSolution,
#                                      move,
#                                      network,
#                                      demand,
#                                      tp;
#                                      full_threshold::Float64 = 0.70)

#     # FULL evaluation for local-search moves on loop lines
#     if move isa Union{
#             VNSInsertionTerminal,
#             VNSRemovalTerminal,
#             VNSInsertion,
#             VNSRemoval,
#             VNSSwapWithinLine,
#             VNSSubstitution,
#             VNSReversal
#         }

#         if is_loop_line(current_sol.lineplan, move.line)
#             sol = build_solution_feasible(lineplan, network, demand, tp)
#             return sol, false, 0.0
#         end
#     end

#     # FULL evaluation for disruptive shaking moves
#     if move isa VNSPartialSegmentSwap
#         sol = build_solution_feasible(lineplan, network, demand, tp)
#         return sol, false, 0.0
#     end

#     # Dijkstra PQ for least disruptive moves
#     if move isa VNSInsertionTerminal ||
#        move isa VNSRemovalTerminal ||
#        move isa VNSReversal

#         sol, detection_time =
#             build_solution_local_dijkstra_pq(
#                 lineplan,
#                 current_sol,
#                 move,
#                 network,
#                 demand,
#                 tp
#             )

#         return sol, false, detection_time
#     end

#     # Threshold basic for medium-disruption moves
#     sol, detection_time, method_used =
#         build_solution_local_threshold_basic(
#             lineplan,
#             current_sol,
#             move,
#             network,
#             demand,
#             tp;
#             full_threshold = full_threshold
#         )

#     return sol, false, detection_time
# end

# UNIFIED SOLUTION BUILDERS

# Build a TransitSolution using the selected evaluation mode
function build_solution(lineplan,
                        current_sol,
                        move,
                        eval_mode::Symbol,
                        network,
                        demand,
                        tp;
                        iteration::Int = 0,
                        full_eval_every::Int = 50)

    if eval_mode == :full
        return build_solution_feasible(lineplan, network, demand, tp)

    elseif eval_mode == :local_dijkstra_basic
        sol, _ = build_solution_local_dijkstra_basic(lineplan, current_sol, move, network, demand, tp)
        return sol

    elseif eval_mode == :local_dijkstra_pq
        sol, _ = build_solution_local_dijkstra_pq(lineplan, current_sol, move, network, demand, tp)
        return sol

    elseif eval_mode == :local_floyd_subset
        sol, _ = build_solution_local_floyd_subset(lineplan, current_sol, move, network, demand, tp)
        return sol

    elseif eval_mode == :local_floyd_subset_path
        sol, _ = build_solution_local_floyd_subset_path(lineplan, current_sol, move, network, demand, tp)
        return sol

    elseif eval_mode == :local_floyd_subset_path_fallback
        sol, _ = build_solution_local_floyd_subset_path_fallback(
            lineplan, current_sol, move, network, demand, tp;
            dijkstra_mode = :pq
        )
        return sol

    elseif eval_mode == :local_fw_subset_path_basic
        sol, _ = build_solution_local_fw_subset_path_basic(lineplan, current_sol, move, network, demand, tp)
        return sol

    elseif eval_mode == :local_fw_subset_path_pq
        sol, _ = build_solution_local_fw_subset_path_pq(lineplan, current_sol, move, network, demand, tp)
        return sol

    elseif eval_mode == :local_hybrid_basic
        sol, _, _ = build_solution_local_hybrid_basic(lineplan, current_sol, move, network, demand, tp)
        return sol

    elseif eval_mode == :local_hybrid_pq
        sol, _, _ = build_solution_local_hybrid_pq(lineplan, current_sol, move, network, demand, tp)
        return sol

    elseif eval_mode == :local_threshold_adaptive
        sol, _, _ = build_solution_local_threshold_adaptive(
            lineplan, current_sol, move, network, demand, tp;
            full_threshold = 0.75,
            dijkstra_mode = :pq
        )
        return sol

    elseif eval_mode == :local_threshold_basic
        sol, _, _ = build_solution_local_threshold_basic(
            lineplan, current_sol, move, network, demand, tp;
            full_threshold = 0.75
        )
        return sol

    elseif eval_mode == :local_threshold_pq
        sol, _, _ = build_solution_local_threshold_pq(
            lineplan, current_sol, move, network, demand, tp;
            full_threshold = 0.75
        )
        return sol
    
    # For adaptive method testing
    # elseif eval_mode == :adaptive_local
    #     sol, _, _ =
    #         build_solution_vns_adaptive(
    #             lineplan,
    #             current_sol,
    #             move,
    #             network,
    #             demand,
    #             tp;
    #             full_threshold = 0.70
    #         )
    #     return sol
        
    elseif eval_mode == :local_floyd_subset_periodic
        if should_do_periodic_full_eval(iteration, full_eval_every)
            return build_solution_feasible(lineplan, network, demand, tp)
        else
            sol, _ = build_solution_local_floyd_subset(lineplan, current_sol, move, network, demand, tp)
            return sol
        end

    elseif eval_mode == :local_hybrid_basic_periodic
        if should_do_periodic_full_eval(iteration, full_eval_every)
            return build_solution_feasible(lineplan, network, demand, tp)
        else
            sol, _, _ = build_solution_local_hybrid_basic(lineplan, current_sol, move, network, demand, tp)
            return sol
        end

    elseif eval_mode == :local_hybrid_pq_periodic
        if should_do_periodic_full_eval(iteration, full_eval_every)
            return build_solution_feasible(lineplan, network, demand, tp)
        else
            sol, _, _ = build_solution_local_hybrid_pq(lineplan, current_sol, move, network, demand, tp)
            return sol
        end

    else
        error("Unknown eval_mode: $eval_mode")
    end
end

# Unified dispatcher for building a TransitSolution with evaluation metadata
function build_solution_with_meta(lineplan,
                                  current_sol,
                                  move,
                                  eval_mode::Symbol,
                                  network,
                                  demand,
                                  tp;
                                  iteration::Int = 0,
                                  full_eval_every::Int = 50)

    if eval_mode == :full
        sol = build_solution_feasible(lineplan, network, demand, tp)
        return sol, false, 0.0

    elseif eval_mode == :local_dijkstra_basic
        sol, detection_time = build_solution_local_dijkstra_basic(lineplan, current_sol, move, network, demand, tp)
        return sol, false, detection_time

    elseif eval_mode == :local_dijkstra_pq
        sol, detection_time = build_solution_local_dijkstra_pq(lineplan, current_sol, move, network, demand, tp)
        return sol, false, detection_time

    elseif eval_mode == :local_floyd_subset
        sol, detection_time = build_solution_local_floyd_subset(lineplan, current_sol, move, network, demand, tp)
        return sol, false, detection_time

    elseif eval_mode == :local_floyd_subset_path
        sol, detection_time = build_solution_local_floyd_subset_path(lineplan, current_sol, move, network, demand, tp)
        return sol, false, detection_time

    elseif eval_mode == :local_floyd_subset_path_fallback
        sol, detection_time = build_solution_local_floyd_subset_path_fallback(
            lineplan, current_sol, move, network, demand, tp;
            dijkstra_mode = :pq
        )
        return sol, false, detection_time

    elseif eval_mode == :local_fw_subset_path_basic
        sol, detection_time = build_solution_local_fw_subset_path_basic(lineplan, current_sol, move, network, demand, tp)
        return sol, false, detection_time

    elseif eval_mode == :local_fw_subset_path_pq
        sol, detection_time = build_solution_local_fw_subset_path_pq(lineplan, current_sol, move, network, demand, tp)
        return sol, false, detection_time

    elseif eval_mode == :local_hybrid_basic
        return build_solution_local_hybrid_basic(lineplan, current_sol, move, network, demand, tp)

    elseif eval_mode == :local_hybrid_pq
        return build_solution_local_hybrid_pq(lineplan, current_sol, move, network, demand, tp)

    elseif eval_mode == :local_threshold_adaptive
        sol, detection_time, _ = build_solution_local_threshold_adaptive(
            lineplan, current_sol, move, network, demand, tp;
            full_threshold = 0.70,
            dijkstra_mode = :pq
        )
        return sol, false, detection_time

    elseif eval_mode == :local_threshold_basic
        sol, detection_time, _ = build_solution_local_threshold_basic(
            lineplan, current_sol, move, network, demand, tp;
            full_threshold = 0.70
        )
        return sol, false, detection_time

    elseif eval_mode == :local_threshold_pq
        sol, detection_time, _ = build_solution_local_threshold_pq(
            lineplan, current_sol, move, network, demand, tp;
            full_threshold = 0.70
        )
        return sol, false, detection_time

    # For adaptive method testing
    # elseif eval_mode == :adaptive_local
    #     return build_solution_vns_adaptive(
    #         lineplan,
    #         current_sol,
    #         move,
    #         network,
    #         demand,
    #         tp;
    #         full_threshold = 0.70
    #     )

    elseif eval_mode == :local_floyd_subset_periodic
        if should_do_periodic_full_eval(iteration, full_eval_every)
            sol = build_solution_feasible(lineplan, network, demand, tp)
            return sol, false, 0.0
        else
            sol, detection_time = build_solution_local_floyd_subset(lineplan, current_sol, move, network, demand, tp)
            return sol, false, detection_time
        end

    elseif eval_mode == :local_hybrid_basic_periodic
        if should_do_periodic_full_eval(iteration, full_eval_every)
            sol = build_solution_feasible(lineplan, network, demand, tp)
            return sol, false, 0.0
        else
            return build_solution_local_hybrid_basic(lineplan, current_sol, move, network, demand, tp)
        end

    elseif eval_mode == :local_hybrid_pq_periodic
        if should_do_periodic_full_eval(iteration, full_eval_every)
            sol = build_solution_feasible(lineplan, network, demand, tp)
            return sol, false, 0.0
        else
            return build_solution_local_hybrid_pq(lineplan, current_sol, move, network, demand, tp)
        end

    else
        error("Unknown eval_mode: $eval_mode")
    end
end