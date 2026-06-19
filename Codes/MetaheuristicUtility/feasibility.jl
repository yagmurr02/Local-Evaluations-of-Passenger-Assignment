include(joinpath(@__DIR__, "common_lineplan_utils.jl"))

# HELPER FUNCTIONS

# Return the set of nodes covered by the lineplan
function covered_nodes(lineplan)
    covered = Set{Int}()
    for L in 1:size(lineplan, 1)
        route = route_from_lineplan(lineplan, L)
        for v in route
            push!(covered, v)
        end
    end
    return covered
end

# Return a Boolean vector indicating which nodes appear in the lineplan
function nodes_present(lineplan, n_nodes::Int)
    present = falses(n_nodes)
    for L in 1:size(lineplan, 1)
        for v in route_from_lineplan(lineplan, L)
            present[v] = true
        end
    end
    return present
end

# Return the set of nodes appearing in positive-demand OD pairs
function demand_nodes(demand, n_nodes::Int)
    nodes = Set{Int}()
    for i in 1:n_nodes
        for j in 1:n_nodes
            if i != j && demand[i, j] > 0
                push!(nodes, i)
                push!(nodes, j)
            end
        end
    end
    return nodes
end

# Return a Boolean vector indicating nodes involved in positive-demand OD pairs
function demand_nodes_present(demand, n_nodes::Int)
    relevant = falses(n_nodes)
    for i in 1:n_nodes
        for j in 1:n_nodes
            if i != j && demand[i, j] > 0
                relevant[i] = true
                relevant[j] = true
            end
        end
    end
    return relevant
end

# Construct the route-set graph induced by the lineplan
function build_route_set_graph(lineplan, n_nodes::Int)
    adj = [Int[] for _ in 1:n_nodes]

    for L in 1:size(lineplan, 1)
        route = route_from_lineplan(lineplan, L)

        if length(route) >= 2
            for i in 1:(length(route) - 1)
                u = route[i]
                v = route[i + 1]
                push!(adj[u], v)
                push!(adj[v], u)
            end
        end
    end

    for i in 1:n_nodes
        adj[i] = unique(adj[i])
    end

    return adj
end

# FEASIBILITY FUNCTIONS

# Check whether the route-set graph is connected over demand-relevant nodes
function is_connected_route_set(lineplan, demand, n_nodes::Int)
    relevant = demand_nodes(demand, n_nodes)
    covered = covered_nodes(lineplan)

    isempty(relevant) && return false

    for v in relevant
        if !(v in covered)
            return false
        end
    end

    g = build_route_set_graph(lineplan, n_nodes)
    start = first(relevant)

    visited = Set{Int}()
    stack = [start]

    while !isempty(stack)
        u = pop!(stack)

        if u in visited
            continue
        end

        push!(visited, u)

        for v in g[u]
            if (v in relevant) && !(v in visited)
                push!(stack, v)
            end
        end
    end

    return length(visited) == length(relevant)
end

# Check whether all positive-demand OD pairs are served by the lineplan
function serves_all_od_pairs(lineplan,
                             network,
                             demand,
                             n_nodes::Int,
                             tp::Real)

    shortest, _, _, _ = shortest_paths(network, lineplan, tp)

    for o in 1:n_nodes
        for d in 1:n_nodes
            if o != d && demand[o, d] > 0 && !isfinite(shortest[o, d])
                return false
            end
        end
    end

    return true
end

# Check whether a lineplan satisfies route, coverage, and connectivity conditions
function is_feasible(lineplan,
                     network,
                     demand,
                     n_nodes::Int;
                     tp::Real = 5.0,
                     min_nodes::Int = 2,
                     max_nodes::Int = typemax(Int),
                     enforce_connected::Bool = false,
                     enforce_all_covered::Bool = true,
                     allow_loop_lines::Bool = true,
                     verbose::Bool = false)

    for L in 1:size(lineplan, 1)
        route = route_from_lineplan(lineplan, L)

        if isempty(route)
            verbose && println("Line $L infeasible: empty route")
            return false
        end

        if !is_valid_route(route, network;
                           min_nodes = min_nodes,
                           max_nodes = max_nodes,
                           allow_loop_lines = allow_loop_lines)
            verbose && println("Line $L infeasible: invalid route structure or missing edge")
            return false
        end
    end

    if enforce_all_covered
        if !serves_all_od_pairs(lineplan, network, demand, n_nodes, tp)
            verbose && println("Lineplan infeasible: not all positive-demand OD-pairs are served")
            return false
        end
    end

    if enforce_connected
        if !is_connected_route_set(lineplan, demand, n_nodes)
            verbose && println("Lineplan infeasible: route-set graph is disconnected over demand-relevant nodes")
            return false
        end
    end

    return true
end