using Random

include(joinpath(@__DIR__, "common_lineplan_utils.jl"))
include(joinpath(@__DIR__, "..", "LocalMoves_BEM.jl"))
include(joinpath(@__DIR__, "basic_vns_moves.jl"))

# Extract the stop sequence
function proposed_line_stops(lineplan, L::Int)
    return route_from_lineplan(lineplan, L)
end

# Get candidate nodes for insertion that are not already in the route
function candidate_nodes_not_in_route(route::Vector{Int}, n_nodes::Int)
    used = Set(route)
    return [v for v in 1:n_nodes if !(v in used)]
end

# Generate all unordered pairs of nodes in the route
function unordered_node_pairs(route::Vector{Int}; require_internal_segment::Bool = false)
    pairs = Tuple{Int,Int}[]
    n = length(route)

    n < 2 && return pairs

    for i in 1:(n - 1)
        for j in (i + 1):n
            if require_internal_segment
                j - i <= 1 && continue
            end
            push!(pairs, (route[i], route[j]))
        end
    end

    return pairs
end

# Generate consecutive pairs of nodes in the route
function consecutive_pairs(route::Vector{Int})
    length(route) < 2 && return Tuple{Int,Int}[]
    return [(route[i], route[i + 1]) for i in 1:(length(route) - 1)]
end

# Check if the candidate lineplan is different from the original
function changed_lineplan(original, candidate)
    candidate === nothing && return false
    return !(candidate == original)
end

# Check if the candidate lineplan is feasible
function candidate_is_feasible(cand,
                               network,
                               demand,
                               tp,
                               n_nodes;
                               min_nodes::Int,
                               max_nodes::Int)

    cand === nothing && return false

    return is_feasible(
        cand,
        network,
        demand,
        n_nodes;
        tp = tp,
        min_nodes = min_nodes,
        max_nodes = max_nodes,
        enforce_connected = false,
        enforce_all_covered = true
    )
end

# Safe version of terminal insertion that return `nothing` if the move is not possible
function safe_insertion_at_terminal(lineplan, L, x; position::String)
    try
        return insertion_at_terminal(lineplan, L, x; position = position)
    catch
        return nothing
    end
end

# Safe version of terminal removal that return `nothing` if the move is not possible
function safe_removal_at_terminal(lineplan, L; position::String)
    try
        return removal_at_terminal(lineplan, L; position = position)
    catch
        return nothing
    end
end

# Insertion at terminal
function proposed_random_terminal_insertion_move(lineplan,
                                                 network,
                                                 rng;
                                                 min_nodes = 2,
                                                 max_nodes = typemax(Int))

    n_routes = size(lineplan, 1)
    n_nodes  = size(network, 1)

    for L in randperm(rng, n_routes)
        route = proposed_line_stops(lineplan, L)

        isempty(route) && continue
        length(route) >= max_nodes && continue

        candidates = candidate_nodes_not_in_route(route, n_nodes)
        isempty(candidates) && continue
        shuffle!(rng, candidates)

        positions = shuffle!(rng, ["start", "end"])

        for pos in positions
            for x in candidates
                cand = safe_insertion_at_terminal(lineplan, L, x; position = pos)
                cand === nothing && continue

                if changed_lineplan(lineplan, cand)
                    move = ProposedInsertionTerminal(L, x, pos)
                    return cand, move, "InsertionTerminal"
                end
            end
        end
    end

    return nothing, nothing, "InsertionTerminal"
end

# Removal at terminal
function proposed_random_terminal_removal_move(lineplan,
                                               network,
                                               rng;
                                               min_nodes = 2,
                                               max_nodes = typemax(Int))

    for L in randperm(rng, size(lineplan, 1))
        route = proposed_line_stops(lineplan, L)

        length(route) <= min_nodes && continue
        length(route) <= 2 && continue

        positions = shuffle!(rng, ["start", "end"])

        for pos in positions
            cand = safe_removal_at_terminal(lineplan, L; position = pos)
            cand === nothing && continue

            if changed_lineplan(lineplan, cand)
                removed_node = pos == "start" ? route[1] : route[end]
                move = ProposedRemovalTerminal(L, removed_node, pos)
                return cand, move, "RemovalTerminal"
            end
        end
    end

    return nothing, nothing, "RemovalTerminal"
end

# Reversal
function proposed_random_reversal_move(lineplan,
                                       network,
                                       rng;
                                       min_nodes = 2,
                                       max_nodes = typemax(Int))

    for L in randperm(rng, size(lineplan, 1))
        route = proposed_line_stops(lineplan, L)

        length(route) < 2 && continue

        cand = reversal(lineplan, L)

        if changed_lineplan(lineplan, cand)
            move = ProposedReversal(L)
            return cand, move, "Reversal"
        end
    end

    return nothing, nothing, "Reversal"
end

# Insertion
function proposed_random_insertion_move(lineplan,
                                        network,
                                        rng;
                                        min_nodes = 2,
                                        max_nodes = typemax(Int))

    n_routes = size(lineplan, 1)
    n_nodes  = size(network, 1)

    for L in randperm(rng, n_routes)
        route = proposed_line_stops(lineplan, L)

        length(route) < 2 && continue
        length(route) >= max_nodes && continue

        segs = consecutive_pairs(route)
        isempty(segs) && continue
        shuffle!(rng, segs)

        candidates = candidate_nodes_not_in_route(route, n_nodes)
        isempty(candidates) && continue
        shuffle!(rng, candidates)

        for (a, b) in segs
            for x in candidates
                cand = insertion(lineplan, L, a, b, x, network)

                if changed_lineplan(lineplan, cand)
                    move = ProposedInsertion(L, a, b, x)
                    return cand, move, "InsertionInternal"
                end
            end
        end
    end

    return nothing, nothing, "InsertionInternal"
end

# Removal
function proposed_random_removal_move(lineplan,
                                      network,
                                      rng;
                                      min_nodes = 2,
                                      max_nodes = typemax(Int))

    for L in randperm(rng, size(lineplan, 1))
        route = proposed_line_stops(lineplan, L)

        length(route) <= max(min_nodes, 2) && continue
        length(route) <= 2 && continue

        internal_nodes = route[2:end-1]
        isempty(internal_nodes) && continue

        for x in shuffle!(rng, copy(internal_nodes))
            cand = removal(lineplan, L, x)

            if changed_lineplan(lineplan, cand)
                move = ProposedRemoval(L, x)
                return cand, move, "RemovalInternal"
            end
        end
    end

    return nothing, nothing, "RemovalInternal"
end

# Swap within line
function proposed_random_swap_within_line_move(lineplan,
                                               network,
                                               rng;
                                               min_nodes = 2,
                                               max_nodes = typemax(Int))

    for L in randperm(rng, size(lineplan, 1))
        route = proposed_line_stops(lineplan, L)

        length(route) < 2 && continue

        pairs = unordered_node_pairs(route)
        isempty(pairs) && continue
        shuffle!(rng, pairs)

        for (a, b) in pairs
            cand = swap_within_line(lineplan, L, a, b, network)

            if changed_lineplan(lineplan, cand)
                move = ProposedSwapWithinLine(L, a, b)
                return cand, move, "SwapWithinLine"
            end
        end
    end

    return nothing, nothing, "SwapWithinLine"
end

# Substitution
function proposed_random_substitution_move(lineplan,
                                           network,
                                           rng;
                                           min_nodes = 2,
                                           max_nodes = typemax(Int))

    n_nodes = size(network, 1)

    for L in randperm(rng, size(lineplan, 1))
        route = proposed_line_stops(lineplan, L)
        isempty(route) && continue

        for old in shuffle!(rng, copy(route))
            candidates = [v for v in candidate_nodes_not_in_route(route, n_nodes) if v != old]
            isempty(candidates) && continue
            shuffle!(rng, candidates)

            for new in candidates
                cand = substitution(lineplan, L, old, new, network)

                if changed_lineplan(lineplan, cand)
                    move = ProposedSubstitution(L, old, new)
                    return cand, move, "Substitution"
                end
            end
        end
    end

    return nothing, nothing, "Substitution"
end

# Segment reversal
function proposed_random_segment_reversal_move(lineplan,
                                               network,
                                               rng;
                                               min_nodes = 2,
                                               max_nodes = typemax(Int))

    for L in randperm(rng, size(lineplan, 1))
        route = proposed_line_stops(lineplan, L)

        length(route) < 3 && continue

        pairs = unordered_node_pairs(route; require_internal_segment = true)
        isempty(pairs) && continue
        shuffle!(rng, pairs)

        for (start_node, end_node) in pairs
            cand = segment_reversal(lineplan, L, start_node, end_node)

            if changed_lineplan(lineplan, cand)
                move = ProposedSegmentReversal(L, start_node, end_node)
                return cand, move, "SegmentReversal"
            end
        end
    end

    return nothing, nothing, "SegmentReversal"
end

# Swap between lines
function proposed_random_swap_between_lines_move(lineplan,
                                                 network,
                                                 rng;
                                                 min_nodes = 2,
                                                 max_nodes = typemax(Int))

    n_routes = size(lineplan, 1)
    n_routes < 2 && return nothing, nothing, "SwapBetweenLines"

    line_pairs = [(i, j) for i in 1:(n_routes - 1) for j in (i + 1):n_routes]
    shuffle!(rng, line_pairs)

    for (La, Lb) in line_pairs
        ra = proposed_line_stops(lineplan, La)
        rb = proposed_line_stops(lineplan, Lb)

        isempty(ra) && continue
        isempty(rb) && continue

        nodes_a = shuffle!(rng, copy(ra))
        nodes_b = shuffle!(rng, copy(rb))

        for a in nodes_a
            for b in nodes_b
                cand = swap_between_lines(lineplan, La, Lb, a, b, network)

                if changed_lineplan(lineplan, cand)
                    move = ProposedSwapBetweenLines(La, a, Lb, b)
                    return cand, move, "SwapBetweenLines"
                end
            end
        end
    end

    return nothing, nothing, "SwapBetweenLines"
end

# Transfer between lines
function proposed_random_position_change_between_lines_move(lineplan,
                                                            network,
                                                            rng;
                                                            min_nodes = 2,
                                                            max_nodes = typemax(Int))

    n_routes = size(lineplan, 1)
    n_routes < 2 && return nothing, nothing, "PositionChangeBetweenLines"

    line_pairs = [(i, j) for i in 1:n_routes for j in 1:n_routes if i != j]
    shuffle!(rng, line_pairs)

    for (from_line, to_line) in line_pairs
        from_route = proposed_line_stops(lineplan, from_line)
        to_route   = proposed_line_stops(lineplan, to_line)

        length(from_route) <= min_nodes && continue
        length(to_route) < 2 && continue
        length(to_route) >= max_nodes && continue

        moved_indices = collect(1:length(from_route))
        shuffle!(rng, moved_indices)

        segs = consecutive_pairs(to_route)
        isempty(segs) && continue
        shuffle!(rng, segs)

        for idx in moved_indices
            moved_node = from_route[idx]

            for (left_node, right_node) in segs
                cand = transfer_between_lines(
                    lineplan,
                    from_line,
                    to_line,
                    moved_node,
                    left_node,
                    right_node,
                    network
                )

                if changed_lineplan(lineplan, cand)
                    move = ProposedPositionChangeBetweenLines(
                        from_line,
                        to_line,
                        moved_node,
                        left_node,
                        right_node
                    )
                    return cand, move, "PositionChangeBetweenLines"
                end
            end
        end
    end

    return nothing, nothing, "PositionChangeBetweenLines"
end

# Partial line swap
function proposed_random_partial_line_swap_move(lineplan,
                                                network,
                                                rng;
                                                min_nodes = 2,
                                                max_nodes = typemax(Int))

    n_routes = size(lineplan, 1)
    n_routes < 2 && return nothing, nothing, "PartialLineSwap"

    line_pairs = [(i, j) for i in 1:(n_routes - 1) for j in (i + 1):n_routes]
    shuffle!(rng, line_pairs)

    for (La, Lb) in line_pairs
        ra = proposed_line_stops(lineplan, La)
        rb = proposed_line_stops(lineplan, Lb)

        length(ra) < 3 && continue
        length(rb) < 3 && continue

        ea = unordered_node_pairs(ra; require_internal_segment = true)
        eb = unordered_node_pairs(rb; require_internal_segment = true)

        isempty(ea) && continue
        isempty(eb) && continue

        shuffle!(rng, ea)
        shuffle!(rng, eb)

        for (aj, ak) in ea
            for (bl, bm) in eb
                cand = partial_segment_swap(lineplan, La, aj, ak, Lb, bl, bm, network)

                if changed_lineplan(lineplan, cand)
                    move = ProposedPartialLineSwap(La, aj, ak, Lb, bl, bm)
                    return cand, move, "PartialLineSwap"
                end
            end
        end
    end

    return nothing, nothing, "PartialLineSwap"
end

# Move lists
const PROPOSED_INTENSIFICATION_FNS = [
    proposed_random_terminal_insertion_move,
    proposed_random_terminal_removal_move,
    proposed_random_insertion_move,
    proposed_random_removal_move,
    proposed_random_swap_within_line_move,
    proposed_random_substitution_move,
    proposed_random_reversal_move,
]

const PROPOSED_SHAKING_NEIGHBORHOODS = [
    proposed_random_swap_between_lines_move,
    proposed_random_position_change_between_lines_move,
    proposed_random_segment_reversal_move,
    proposed_random_partial_line_swap_move,
]

# Intensification neighbor
function proposed_random_intensification_neighbor(lineplan,
                                                  network,
                                                  demand,
                                                  tp,
                                                  n_nodes::Int,
                                                  rng;
                                                  min_nodes = 2,
                                                  max_nodes = typemax(Int))

    for fn in PROPOSED_INTENSIFICATION_FNS

        cand, move, name = fn(
            lineplan,
            network,
            rng;
            min_nodes = min_nodes,
            max_nodes = max_nodes
        )

        cand === nothing && continue
        changed_lineplan(lineplan, cand) || continue

        if candidate_is_feasible(
                cand,
                network,
                demand,
                tp,
                n_nodes;
                min_nodes = min_nodes,
                max_nodes = max_nodes
           )
            return cand, move, name
        end
    end

    return nothing, nothing, "None"
end

# Shaking count
function proposed_shaking_neighbourhood_count()
    return length(PROPOSED_SHAKING_NEIGHBORHOODS)
end

# Shaking neighbor
function proposed_random_shake_neighbor(lineplan,
                                        k::Int,
                                        network,
                                        demand,
                                        tp,
                                        n_nodes::Int,
                                        rng;
                                        min_nodes = 2,
                                        max_nodes = typemax(Int))

    if k < 1 || k > proposed_shaking_neighbourhood_count()
        return nothing, nothing, "None"
    end

    fn = PROPOSED_SHAKING_NEIGHBORHOODS[k]

    cand, move, name = fn(
        lineplan,
        network,
        rng;
        min_nodes = min_nodes,
        max_nodes = max_nodes
    )

    cand === nothing && return nothing, nothing, name
    changed_lineplan(lineplan, cand) || return nothing, nothing, name

    if candidate_is_feasible(
            cand,
            network,
            demand,
            tp,
            n_nodes;
            min_nodes = min_nodes,
            max_nodes = max_nodes
       )
        return cand, move, name
    end

    return nothing, nothing, name
end

# Generate feasible proposed intensification neighbor with retries
function generate_feasible_proposed_intensification_neighbor(sol::TransitSolution,
                                                             network,
                                                             demand,
                                                             tp,
                                                             n_nodes::Int,
                                                             rng::AbstractRNG;
                                                             min_nodes::Int = 2,
                                                             max_nodes::Int = typemax(Int),
                                                             max_tries::Int = 10_000)

    for _ in 1:max_tries
        cand, move, name = proposed_random_intensification_neighbor(
            sol.lineplan,
            network,
            demand,
            tp,
            n_nodes,
            rng;
            min_nodes = min_nodes,
            max_nodes = max_nodes
        )

        cand === nothing && continue
        changed_lineplan(sol.lineplan, cand) || continue

        return cand, move, name
    end

    return nothing, nothing, "None"
end

# Generate feasible proposed shake neighbor with retries
function generate_feasible_proposed_shake_neighbor(sol::TransitSolution,
                                                   k::Int,
                                                   network,
                                                   demand,
                                                   tp,
                                                   n_nodes::Int,
                                                   rng::AbstractRNG;
                                                   min_nodes::Int = 2,
                                                   max_nodes::Int = typemax(Int),
                                                   max_tries::Int = 20)

    for _ in 1:max_tries
        cand, move, name = proposed_random_shake_neighbor(
            sol.lineplan,
            k,
            network,
            demand,
            tp,
            n_nodes,
            rng;
            min_nodes = min_nodes,
            max_nodes = max_nodes
        )

        cand === nothing && continue
        changed_lineplan(sol.lineplan, cand) || continue

        return cand, move, name
    end

    return nothing, nothing, "None"
end