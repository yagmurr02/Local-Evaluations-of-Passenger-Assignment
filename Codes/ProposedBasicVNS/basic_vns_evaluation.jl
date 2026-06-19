include(joinpath(@__DIR__, "basic_vns_moves.jl"))
include(joinpath(@__DIR__, "evaluation.jl"))

# Detect affected OD pairs for insertion at terminal moves
function detect_affected_ods(current_sol::TransitSolution,
                             cand_lineplan,
                             move::ProposedInsertionTerminal,
                             network,
                             tp)

    directlink_new, lineinfo_new = directlinknetworkinit(network, cand_lineplan)

    return detect_insertion_terminal(
        current_sol.shortest,
        current_sol.directlink,
        directlink_new,
        current_sol.lineplan,
        move.line,
        move.node,
        tp;
        path_old = current_sol.path_old,
        lineinfo_old = current_sol.lineinfo,
        lineinfo_new = lineinfo_new
    )
end

# Detect affected OD pairs for removal at terminal moves
function detect_affected_ods(current_sol::TransitSolution,
                             cand_lineplan,
                             move::ProposedRemovalTerminal,
                             network,
                             tp)

    directlink_new, lineinfo_new = directlinknetworkinit(network, cand_lineplan)

    return detect_removal_terminal(
        current_sol.shortest,
        current_sol.directlink,
        directlink_new,
        current_sol.lineplan,
        move.line,
        move.node,
        tp;
        path_old = current_sol.path_old,
        lineinfo_old = current_sol.lineinfo,
        lineinfo_new = lineinfo_new
    )
end

# Detect affected OD pairs for reversal moves
function detect_affected_ods(current_sol::TransitSolution,
                             cand_lineplan,
                             move::ProposedReversal,
                             network,
                             tp)

    directlink_new, lineinfo_new = directlinknetworkinit(network, cand_lineplan)

    return detect_reversal(
        current_sol.shortest,
        current_sol.directlink,
        directlink_new,
        current_sol.lineplan,
        move.line,
        tp;
        path_old = current_sol.path_old,
        lineinfo_old = current_sol.lineinfo,
        lineinfo_new = lineinfo_new
    )
end

# Detect affected OD pairs for insertion moves
function detect_affected_ods(current_sol::TransitSolution,
                             cand_lineplan,
                             move::ProposedInsertion,
                             network,
                             tp)

    directlink_new, lineinfo_new = directlinknetworkinit(network, cand_lineplan)

    return detect_insertion_internal(
        current_sol.shortest,
        current_sol.directlink,
        directlink_new,
        current_sol.lineplan,
        move.line,
        move.left_node,
        move.right_node,
        move.inserted_node,
        tp;
        path_old = current_sol.path_old,
        lineinfo_old = current_sol.lineinfo,
        lineinfo_new = lineinfo_new
    )
end

# Detect affected OD pairs for removal moves
function detect_affected_ods(current_sol::TransitSolution,
                             cand_lineplan,
                             move::ProposedRemoval,
                             network,
                             tp)

    directlink_new, lineinfo_new = directlinknetworkinit(network, cand_lineplan)

    return detect_removal_internal(
        current_sol.shortest,
        current_sol.directlink,
        directlink_new,
        current_sol.lineplan,
        move.line,
        move.removed_node,
        tp;
        path_old = current_sol.path_old,
        lineinfo_old = current_sol.lineinfo,
        lineinfo_new = lineinfo_new
    )
end

# Detect affected OD pairs for swap within line moves
function detect_affected_ods(current_sol::TransitSolution,
                             cand_lineplan,
                             move::ProposedSwapWithinLine,
                             network,
                             tp)

    directlink_new, lineinfo_new = directlinknetworkinit(network, cand_lineplan)

    return detect_swap_within_line(
        current_sol.shortest,
        current_sol.directlink,
        directlink_new,
        current_sol.lineplan,
        move.line,
        move.node1,
        move.node2,
        tp;
        path_old = current_sol.path_old,
        lineinfo_old = current_sol.lineinfo,
        lineinfo_new = lineinfo_new
    )
end

# Detect affected OD pairs for substitution moves
function detect_affected_ods(current_sol::TransitSolution,
                             cand_lineplan,
                             move::ProposedSubstitution,
                             network,
                             tp)

    directlink_new, lineinfo_new = directlinknetworkinit(network, cand_lineplan)

    return detect_substitution(
        current_sol.shortest,
        current_sol.directlink,
        directlink_new,
        current_sol.lineplan,
        move.line,
        move.old_node,
        move.new_node,
        tp;
        path_old = current_sol.path_old,
        lineinfo_old = current_sol.lineinfo,
        lineinfo_new = lineinfo_new
    )
end

# Detect affected OD pairs for segment reversal moves
function detect_affected_ods(current_sol::TransitSolution,
                             cand_lineplan,
                             move::ProposedSegmentReversal,
                             network,
                             tp)

    directlink_new, lineinfo_new = directlinknetworkinit(network, cand_lineplan)

    return detect_segment_reversal(
        current_sol.shortest,
        current_sol.directlink,
        directlink_new,
        current_sol.lineplan,
        move.line,
        move.start_node,
        move.end_node,
        tp;
        path_old = current_sol.path_old,
        lineinfo_old = current_sol.lineinfo,
        lineinfo_new = lineinfo_new
    )
end

# Detect affected OD pairs for swap between lines moves
function detect_affected_ods(current_sol::TransitSolution,
                             cand_lineplan,
                             move::ProposedSwapBetweenLines,
                             network,
                             tp)

    directlink_new, lineinfo_new = directlinknetworkinit(network, cand_lineplan)

    return detect_swap_between_lines(
        current_sol.shortest,
        current_sol.directlink,
        directlink_new,
        current_sol.lineplan,
        move.line1,
        move.line2,
        move.node1,
        move.node2,
        tp;
        path_old = current_sol.path_old,
        lineinfo_old = current_sol.lineinfo,
        lineinfo_new = lineinfo_new
    )
end

# Detect affected OD pairs for transfer between lines moves
function detect_affected_ods(current_sol::TransitSolution,
                             cand_lineplan,
                             move::ProposedPositionChangeBetweenLines,
                             network,
                             tp)

    directlink_new, lineinfo_new = directlinknetworkinit(network, cand_lineplan)

    return detect_transfer_between_lines(
        current_sol.shortest,
        current_sol.directlink,
        directlink_new,
        current_sol.lineplan,
        move.from_line,
        move.to_line,
        move.moved_node,
        move.left_node,
        move.right_node,
        tp;
        path_old = current_sol.path_old,
        lineinfo_old = current_sol.lineinfo,
        lineinfo_new = lineinfo_new
    )
end

# Detect affected OD pairs for partial line swap moves
function detect_affected_ods(current_sol::TransitSolution,
                             cand_lineplan,
                             move::ProposedPartialLineSwap,
                             network,
                             tp)

    directlink_new, lineinfo_new = directlinknetworkinit(network, cand_lineplan)

    return detect_partial_segment_swap(
        current_sol.shortest,
        current_sol.directlink,
        directlink_new,
        current_sol.lineplan,
        move.line1,
        move.node1a,
        move.node1b,
        move.line2,
        move.node2a,
        move.node2b,
        tp;
        path_old = current_sol.path_old,
        lineinfo_old = current_sol.lineinfo,
        lineinfo_new = lineinfo_new
    )
end

# Collect stops from the affected line before and after the move for subset construction
function line_stops_for_subset(lineplan_old,
                               cand_lineplan,
                               move::Union{
                                   ProposedInsertion,
                                   ProposedInsertionTerminal,
                                   ProposedRemoval,
                                   ProposedRemovalTerminal,
                                   ProposedSwapWithinLine,
                                   ProposedSubstitution,
                                   ProposedSegmentReversal,
                                   ProposedReversal
                               })

    nodes = Set{Int}()

    old_route = route_from_lineplan(lineplan_old, move.line)
    new_route = route_from_lineplan(cand_lineplan, move.line)

    for v in old_route
        push!(nodes, v)
    end

    for v in new_route
        push!(nodes, v)
    end

    return nodes
end

# Collect stops from both affected lines before and after the move for subset construction
function line_stops_for_subset(lineplan_old,
                               cand_lineplan,
                               move::Union{
                                   ProposedSwapBetweenLines,
                                   ProposedPartialLineSwap
                               })

    nodes = Set{Int}()

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

    return nodes
end

# Collect stops from the source and destination lines before and after the move for subset construction
function line_stops_for_subset(lineplan_old,
                               cand_lineplan,
                               move::ProposedPositionChangeBetweenLines)

    nodes = Set{Int}()

    for line_id in (move.from_line, move.to_line)
        old_route = route_from_lineplan(lineplan_old, line_id)
        new_route = route_from_lineplan(cand_lineplan, line_id)

        for v in old_route
            push!(nodes, v)
        end

        for v in new_route
            push!(nodes, v)
        end
    end

    return nodes
end

# Loop line check helper
function is_loop_line(lineplan, line::Int)
    route = route_from_lineplan(lineplan, line)
    return length(route) != length(unique(route))
end

# Adaptive proposed evaluation selector
function build_solution_proposed_adaptive(lineplan,
                                          current_sol::TransitSolution,
                                          move,
                                          network,
                                          demand,
                                          tp;
                                          full_threshold::Float64 = 0.70)

    if move isa Union{
            ProposedInsertionTerminal,
            ProposedRemovalTerminal,
            ProposedInsertion,
            ProposedRemoval,
            ProposedSwapWithinLine,
            ProposedSubstitution,
            ProposedReversal
        }

        if is_loop_line(current_sol.lineplan, move.line)
            sol = build_solution_feasible(lineplan, network, demand, tp)
            return sol, false, 0.0, :full
        end
    end

    if move isa ProposedSwapBetweenLines ||
        move isa ProposedPositionChangeBetweenLines ||
        move isa ProposedPartialLineSwap ||
        move isa ProposedSegmentReversal

        sol = build_solution_feasible(lineplan, network, demand, tp)
        return sol, false, 0.0, :full
    end

    if move isa ProposedInsertionTerminal ||
        move isa ProposedRemovalTerminal ||
        move isa ProposedReversal

        sol, detection_time =
            build_solution_local_dijkstra_pq(
                lineplan,
                current_sol,
                move,
                network,
                demand,
                tp
            )

        return sol, false, detection_time, :local_dijkstra_pq
    end

    sol, detection_time, method_used =
        build_solution_local_threshold_basic(
            lineplan,
            current_sol,
            move,
            network,
            demand,
            tp;
            full_threshold = full_threshold
        )

    return sol, false, detection_time, method_used
end

# General wrapper for proposed framework
function build_solution_proposed_with_meta(lineplan,
                                           current_sol::TransitSolution,
                                           move,
                                           eval_mode::Symbol,
                                           network,
                                           demand,
                                           tp;
                                           full_threshold::Float64 = 0.70)

    if eval_mode == :full
        sol = build_solution_feasible(lineplan, network, demand, tp)
        return sol, false, 0.0, :full

    elseif eval_mode == :local_dijkstra_pq
        sol, detection_time =
            build_solution_local_dijkstra_pq(
                lineplan,
                current_sol,
                move,
                network,
                demand,
                tp
            )

        return sol, false, detection_time, :local_dijkstra_pq

    elseif eval_mode == :local_floyd_subset
        sol, detection_time =
            build_solution_local_floyd_subset(
                lineplan,
                current_sol,
                move,
                network,
                demand,
                tp
            )

        return sol, false, detection_time, :local_floyd_subset

    elseif eval_mode == :threshold_basic
        sol, detection_time, method_used =
            build_solution_local_threshold_basic(
                lineplan,
                current_sol,
                move,
                network,
                demand,
                tp;
                full_threshold = full_threshold
            )

        return sol, false, detection_time, method_used

    elseif eval_mode == :adaptive_local
        return build_solution_proposed_adaptive(
            lineplan,
            current_sol,
            move,
            network,
            demand,
            tp;
            full_threshold = full_threshold
        )

    else
        error("Unknown proposed evaluation mode: $(eval_mode)")
    end
end