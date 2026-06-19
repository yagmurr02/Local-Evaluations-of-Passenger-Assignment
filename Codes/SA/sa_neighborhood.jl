using Random

include(joinpath(@__DIR__, "common_lineplan_utils.jl"))
include(joinpath(@__DIR__, "..", "LocalMoves_BEM.jl"))

# Terminal insertions: find candidates for insertion at the end of a route
function candidate_terminal_insertions(route::Vector{Int}, network)
    lastnode = route[end]
    used = Set(route)
    candidates = Int[]

    for v in 1:size(network, 1)
        if has_edge(network, lastnode, v) && !(v in used)
            push!(candidates, v)
        end
    end

    return candidates
end

# Check if the candidate lineplan is different from the original
function changed_lineplan(original, candidate)
    return !(candidate == original)
end

# Generate a neighboring lineplan by making a small change to a randomly selected route.
function make_small_change(lineplan,
                           network,
                           rng::AbstractRNG;
                           min_nodes::Int = 2,
                           max_nodes::Int = typemax(Int))

    n_routes = size(lineplan, 1)
    route_order = randperm(rng, n_routes)

    for L in route_order
        route = route_from_lineplan(lineplan, L)
        lenL = length(route)

        if min_nodes < lenL < max_nodes
            if rand(rng) < 0.5
                candidates = candidate_terminal_insertions(route, network)
                if !isempty(candidates)
                    node = rand(rng, candidates)
                    cand = insertion_at_terminal(lineplan, L, node; position = "end")
                    if changed_lineplan(lineplan, cand)
                        return cand, TerminalInsertion(L, node)
                    end
                else
                    cand = reversal(lineplan, L)
                    if changed_lineplan(lineplan, cand)
                        return cand, Reversal(L)
                    end
                end
            else
                cand = removal_at_terminal(lineplan, L; position = "start")
                if changed_lineplan(lineplan, cand)
                    return cand, TerminalRemoval(L)
                end
            end

        elseif lenL == max_nodes
            cand = removal_at_terminal(lineplan, L; position = "start")
            if changed_lineplan(lineplan, cand)
                return cand, TerminalRemoval(L)
            end

        elseif lenL == min_nodes
            candidates = candidate_terminal_insertions(route, network)
            if !isempty(candidates)
                node = rand(rng, candidates)
                cand = insertion_at_terminal(lineplan, L, node; position = "end")
                if changed_lineplan(lineplan, cand)
                    return cand, TerminalInsertion(L, node)
                end
            else
                cand = reversal(lineplan, L)
                if changed_lineplan(lineplan, cand)
                    return cand, Reversal(L)
                end
            end
        end
    end

    return nothing, nothing
end

# Generate a feasible neighboring lineplan by applying small changes to the current solution
function generate_feasible_neighbor(sol::TransitSolution,
                                    network,
                                    demand,
                                    tp,
                                    n_nodes::Int,
                                    rng::AbstractRNG;
                                    min_nodes::Int = 2,
                                    max_nodes::Int = typemax(Int),
                                    max_tries::Int = 10_000)

    for _ in 1:max_tries
        cand_lineplan, move = make_small_change(sol.lineplan, network, rng;
                                                min_nodes = min_nodes,
                                                max_nodes = max_nodes)

        cand_lineplan === nothing && continue

        if is_feasible(cand_lineplan, network, demand, n_nodes;
                       tp = tp,
                       min_nodes = min_nodes,
                       max_nodes = max_nodes,
                       enforce_connected = false,
                       enforce_all_covered = true)
            return cand_lineplan, move
        end
    end

    return nothing, nothing
end 