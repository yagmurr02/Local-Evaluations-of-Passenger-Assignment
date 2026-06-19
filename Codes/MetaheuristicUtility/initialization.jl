include(joinpath(@__DIR__, "common_lineplan_utils.jl"))

using Random

# HELPER FUNCTIONS

# Return all neighboring nodes connected to node u
function neighbors_in_network(network, u::Int)
    n = size(network, 1)
    nbrs = Int[]
    @inbounds for v in 1:n
        if has_edge(network, u, v)
            push!(nbrs, v)
        end
    end
    return nbrs
end

# Identify OD pairs with positive demand but no feasible path
function problematic_od_pairs(lineplan::AbstractMatrix{Int},
                              network,
                              demand,
                              tp)
    n_nodes = size(network, 1)
    shortest, _, _, _ = shortest_paths(network, lineplan, tp)

    bad = Tuple{Int,Int}[]
    for o in 1:n_nodes
        for d in 1:n_nodes
            if o != d && demand[o, d] > 0 && !isfinite(shortest[o, d])
                push!(bad, (o, d))
            end
        end
    end
    return bad
end

# FEASIBILITY CHECKS

# Check whether a lineplan satisfies feasibility conditions
function is_feasible_lineplan(lineplan::AbstractMatrix{Int},
                              network,
                              demand,
                              tp;
                              min_nodes::Int = 2,
                              max_nodes::Int = size(lineplan, 2),
                              allow_loop_lines::Bool = true)

    n_nodes = size(network, 1)

    return is_feasible(lineplan, network, demand, n_nodes;
                       tp = tp,
                       min_nodes = min_nodes,
                       max_nodes = max_nodes,
                       enforce_connected = false,
                       enforce_all_covered = true,
                       allow_loop_lines = allow_loop_lines)
end

# RANDOM ROUTE CONSTRUCTION

# Generate a random connected route satisfying route constraints
function random_connected_route(network,
                                rng::AbstractRNG;
                                min_nodes::Int = 2,
                                max_nodes::Int = 8,
                                preferred_nodes::Union{Nothing,Vector{Int}} = nothing,
                                forced_start::Union{Nothing,Int} = nothing,
                                allow_loop_lines::Bool = true,
                                max_attempts::Int = 200)

    n_nodes = size(network, 1)

    for _ in 1:max_attempts
        target_len = rand(rng, min_nodes:max_nodes)

        start =
            if forced_start !== nothing
                forced_start
            elseif preferred_nodes !== nothing && !isempty(preferred_nodes) && rand(rng) < 0.7
                rand(rng, preferred_nodes)
            else
                rand(rng, 1:n_nodes)
            end

        route = [start]
        used = Set(route)

        while length(route) < target_len
            u = route[end]

            if allow_loop_lines && length(route) == target_len - 1
                cand = [v for v in neighbors_in_network(network, u) if !(v in used) || v == start]
            else
                cand = [v for v in neighbors_in_network(network, u) if !(v in used)]
            end

            isempty(cand) && break

            if preferred_nodes !== nothing
                pref = [v for v in cand if v in preferred_nodes]
                nxt = (!isempty(pref) && rand(rng) < 0.7) ? rand(rng, pref) : rand(rng, cand)
            else
                nxt = rand(rng, cand)
            end

            push!(route, nxt)

            if nxt != start || length(route) < target_len
                push!(used, nxt)
            end
        end

        if is_valid_route(route, network;
                          min_nodes = min_nodes,
                          max_nodes = max_nodes,
                          allow_loop_lines = allow_loop_lines)
            return route
        end
    end

    return nothing
end

# Construct an initial lineplan using randomized route generation
function construct_initial_lineplan(network,
                                    demand,
                                    n_lines::Int,
                                    rng::AbstractRNG;
                                    min_nodes::Int = 2,
                                    max_nodes::Int = 8,
                                    allow_loop_lines::Bool = true,
                                    max_route_attempts::Int = 300)

    n_nodes = size(network, 1)
    relevant = demand_nodes_present(demand, n_nodes)

    lineplan = zeros(Int, n_lines, max_nodes)
    covered = falses(n_nodes)

    relevant_nodes = Int.(findall(relevant))

    first_route = random_connected_route(network, rng;
                                         min_nodes = min_nodes,
                                         max_nodes = max_nodes,
                                         preferred_nodes = isempty(relevant_nodes) ? nothing : relevant_nodes,
                                         forced_start = isempty(relevant_nodes) ? nothing : rand(rng, relevant_nodes),
                                         allow_loop_lines = allow_loop_lines,
                                         max_attempts = max_route_attempts)
    first_route === nothing && return nothing

    set_route!(lineplan, 1, first_route)
    for node in first_route
        covered[node] = true
    end

    for r in 2:n_lines
        uncovered = Int.(findall(relevant .& .!covered))
        covered_nodes_vec = Int.(findall(covered))

        route = nothing

        for _ in 1:max_route_attempts
            forced_start =
                if !isempty(covered_nodes_vec) && rand(rng) < 0.8
                    rand(rng, covered_nodes_vec)
                else
                    nothing
                end

            route = random_connected_route(network, rng;
                                           min_nodes = min_nodes,
                                           max_nodes = max_nodes,
                                           preferred_nodes = isempty(uncovered) ? nothing : uncovered,
                                           forced_start = forced_start,
                                           allow_loop_lines = allow_loop_lines,
                                           max_attempts = 20)
            route !== nothing && break
        end

        route === nothing && return nothing

        set_route!(lineplan, r, route)
        for node in route
            covered[node] = true
        end
    end

    return lineplan
end

# OD COVERAGE REPAIR

# Attempt to repair infeasible OD coverage by replacing a line route
function try_repair_problematic_od(lineplan::AbstractMatrix{Int},
                                   network,
                                   demand,
                                   tp,
                                   rng::AbstractRNG;
                                   min_nodes::Int = 2,
                                   max_nodes::Int = size(lineplan, 2),
                                   allow_loop_lines::Bool = true,
                                   max_attempts::Int = 200)

    bad = problematic_od_pairs(lineplan, network, demand, tp)
    isempty(bad) && return copy(lineplan)

    (o, d) = rand(rng, bad)

    n_nodes = size(network, 1)
    present = nodes_present(lineplan, n_nodes)
    present_nodes_vec = Int.(findall(present))

    preferred_nodes = unique(vcat([o, d], present_nodes_vec))

    for _ in 1:max_attempts
        forced_start =
            if present[o]
                o
            elseif !isempty(present_nodes_vec) && rand(rng) < 0.7
                rand(rng, present_nodes_vec)
            else
                o
            end

        new_route = random_connected_route(network, rng;
                                           min_nodes = min_nodes,
                                           max_nodes = max_nodes,
                                           preferred_nodes = preferred_nodes,
                                           forced_start = forced_start,
                                           allow_loop_lines = allow_loop_lines,
                                           max_attempts = 20)

        new_route === nothing && continue

        candidate = copy(lineplan)
        rid = rand(rng, 1:size(lineplan, 1))
        set_route!(candidate, rid, new_route)

        ok = true
        for r in 1:size(candidate, 1)
            route = line_stops(candidate, r)
            if !is_valid_route(route, network;
                               min_nodes = min_nodes,
                               max_nodes = max_nodes,
                               allow_loop_lines = allow_loop_lines)
                ok = false
                break
            end
        end

        ok || continue

        old_bad = length(bad)
        new_bad = length(problematic_od_pairs(candidate, network, demand, tp))
        if new_bad <= old_bad
            return candidate
        end
    end

    return nothing
end

# Repair OD coverage violations in a lineplan through iterative route replacement
function repair_od_coverage(lineplan::AbstractMatrix{Int},
                            network,
                            demand,
                            tp,
                            rng::AbstractRNG;
                            min_nodes::Int = 2,
                            max_nodes::Int = size(lineplan, 2),
                            allow_loop_lines::Bool = true,
                            max_cover_steps::Int = 200,
                            max_connect_steps::Int = 200)

    current = copy(lineplan)

    for _ in 1:max_cover_steps
        if is_feasible_lineplan(current, network, demand, tp;
                                min_nodes = min_nodes,
                                max_nodes = max_nodes,
                                allow_loop_lines = allow_loop_lines)
            return current
        end

        candidate = try_repair_problematic_od(current, network, demand, tp, rng;
                                              min_nodes = min_nodes,
                                              max_nodes = max_nodes,
                                              allow_loop_lines = allow_loop_lines,
                                              max_attempts = 100)

        candidate === nothing && continue
        current = candidate
    end

    return is_feasible_lineplan(current, network, demand, tp;
                                min_nodes = min_nodes,
                                max_nodes = max_nodes,
                                allow_loop_lines = allow_loop_lines) ? current : nothing
end

# FEASIBLE SOLUTION INITIALIZATION

# Generate a feasible initial lineplan using construction and repair
function initialize_feasible_lineplan(network,
                                      demand,
                                      tp,
                                      n_lines::Int,
                                      rng::AbstractRNG;
                                      min_nodes::Int = 2,
                                      max_nodes::Int = 8,
                                      allow_loop_lines::Bool = true,
                                      outer_attempts::Int = 300,
                                      max_route_attempts::Int = 300,
                                      max_cover_steps::Int = 200,
                                      max_connect_steps::Int = 200,
                                      verbose::Bool = false)

    for attempt in 1:outer_attempts
        lineplan = construct_initial_lineplan(network, demand, n_lines, rng;
                                              min_nodes = min_nodes,
                                              max_nodes = max_nodes,
                                              allow_loop_lines = allow_loop_lines,
                                              max_route_attempts = max_route_attempts)
        lineplan === nothing && continue

        repaired = repair_od_coverage(lineplan, network, demand, tp, rng;
                                      min_nodes = min_nodes,
                                      max_nodes = max_nodes,
                                      allow_loop_lines = allow_loop_lines,
                                      max_cover_steps = max_cover_steps,
                                      max_connect_steps = max_connect_steps)
        repaired === nothing && continue

        if is_feasible_lineplan(repaired, network, demand, tp;
                                min_nodes = min_nodes,
                                max_nodes = max_nodes,
                                allow_loop_lines = allow_loop_lines)
            if verbose
                println("Initialization succeeded on attempt $attempt")
            end
            return repaired
        end
    end

    error("Failed to generate a feasible initial lineplan after $outer_attempts attempts.")
end

# Generate an initial feasible TransitSolution object
function initialize_feasible_solution(network,
                                      demand,
                                      tp,
                                      n_lines::Int,
                                      rng::AbstractRNG;
                                      min_nodes::Int = 2,
                                      max_nodes::Int = 8,
                                      allow_loop_lines::Bool = true,
                                      outer_attempts::Int = 300,
                                      max_route_attempts::Int = 300,
                                      max_cover_steps::Int = 200,
                                      max_connect_steps::Int = 200,
                                      verbose::Bool = false)

    lineplan = initialize_feasible_lineplan(network, demand, tp, n_lines, rng;
                                            min_nodes = min_nodes,
                                            max_nodes = max_nodes,
                                            allow_loop_lines = allow_loop_lines,
                                            outer_attempts = outer_attempts,
                                            max_route_attempts = max_route_attempts,
                                            max_cover_steps = max_cover_steps,
                                            max_connect_steps = max_connect_steps,
                                            verbose = verbose)

    return TransitSolution(lineplan)
end