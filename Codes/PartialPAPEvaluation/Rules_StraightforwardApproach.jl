# Move-specific detection rules for identifying affected OD pairs based on direct-link changes for the straightforward approach.

# HELPERS
# Check whether a line is a loop line
function is_loop_line(stops::Vector{Int})
    return length(stops) >= 2 && first(stops) == last(stops)
end

# For loop lines, get the core stop sequence without the duplicate terminal stop
function loop_core_stops(stops::Vector{Int})
    if is_loop_line(stops)
        return stops[1:end-1]
    else
        return stops
    end
end

# Get stop sequence for a line from lineplan (filtering out zeros)
function line_stops(lineplan, line_id::Int)
    return Int.(filter(!=(0), lineplan[line_id, :]))
end

# Convert changed direct links into endpoint OD pairs under symmetric network assumption
function direct_link_endpoint_pairs_symmetric(edges)
    A = Set{Tuple{Int,Int}}()

    for (u, v) in edges
        u == v && continue
        push!(A, (u, v))
        push!(A, (v, u))
    end

    return A
end

# Convert shortest path matrix to tp-space by adding tp to all finite entries
function shortest_to_tp_space(shortest_old::AbstractMatrix{Float64}, tp::Number)
    n = size(shortest_old, 1)
    shortest_tp = copy(shortest_old)
    @inbounds for i in 1:n, j in 1:n
        if i != j && isfinite(shortest_tp[i,j])
            shortest_tp[i,j] += tp
        end
    end
    return shortest_tp
end

# RULES
# Insertion at terminal
function detect_insertion_terminal(shortest_old::AbstractMatrix{Float64},
                                   directlink_old::AbstractMatrix{Float64},
                                   directlink_new::AbstractMatrix{Float64},
                                   lineplan_old,
                                   line_id::Int,
                                   x::Int,
                                   tp::Number;
                                   tol::Float64 = 1e-9,
                                   path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                   lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                   lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    stops = line_stops(lineplan_old, line_id)
    isempty(stops) && return Set{Tuple{Int,Int}}()

    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    improve_edges = Set{Tuple{Int,Int}}()
    lineonly_changed_edges = Set{Tuple{Int,Int}}()

    @inbounds for s in stops
        # Check incident edge s -> x
        old = directlink_old[s, x]
        new = directlink_new[s, x]

        # Edge detected if newly created or improved
        if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
            push!(improve_edges, (s, x))
        end

        # Edge detected if only the assigned line changes
        if has_lineinfo &&
           isfinite(old) && isfinite(new) &&
           abs(old - new) <= tol &&
           lineinfo_old[s, x] != lineinfo_new[s, x]
            push!(lineonly_changed_edges, (s, x))
        end

        # Check incident edge x → s
        old = directlink_old[x, s]
        new = directlink_new[x, s]

        # Edge detected if newly created or improved
        if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
            push!(improve_edges, (x, s))
        end
        
        # Edge detected if only the assigned line changes
        if has_lineinfo &&
           isfinite(old) && isfinite(new) &&
           abs(old - new) <= tol &&
           lineinfo_old[x, s] != lineinfo_new[x, s]
            push!(lineonly_changed_edges, (x, s))
        end
    end

    isempty(improve_edges) && isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    changed_edges = union(improve_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)
end

# Line addition 
function detect_line_addition(shortest_old::AbstractMatrix{Float64},
                              directlink_old::AbstractMatrix{Float64},
                              directlink_new::AbstractMatrix{Float64},
                              n_vec::Vector{Int},
                              tp::Number;
                              tol::Float64 = 1e-9,
                              path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                              lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                              lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    isempty(n_vec) && return Set{Tuple{Int,Int}}()
    n = size(shortest_old, 1)

    S3 = unique([x for x in n_vec if 1 <= x <= n && x != 0])
    length(S3) < 2 && return Set{Tuple{Int,Int}}()

    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    improve_edges = Set{Tuple{Int,Int}}()

    lineonly_changed_edges = Set{Tuple{Int,Int}}()

    @inbounds for u in S3
        for v in S3
            u == v && continue

            # Check candidate edge u -> v
            old = directlink_old[u, v]
            new = directlink_new[u, v]

            # Edge detected if newly created or improved
            if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
                push!(improve_edges, (u, v))
            end

            # Edge detected if only the assigned line changes
            if has_lineinfo &&
               isfinite(old) && isfinite(new) &&
               abs(old - new) <= tol &&
               lineinfo_old[u, v] != lineinfo_new[u, v]
                push!(lineonly_changed_edges, (u, v))
            end
        end
    end

    isempty(improve_edges) && isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    changed_edges = union(improve_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)
end

# Removal at terminal
function detect_removal_terminal(shortest_old::AbstractMatrix{Float64},
                                 directlink_old::AbstractMatrix{Float64},
                                 directlink_new::AbstractMatrix{Float64},
                                 lineplan_old,
                                 line_id::Int,
                                 y::Int,
                                 tp::Number;
                                 tol::Float64 = 1e-9,
                                 path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                 lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                 lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    stops = line_stops(lineplan_old, line_id)
    isempty(stops) && return Set{Tuple{Int,Int}}()

    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    # Loop-line case
    if is_loop_line(stops)
        core = loop_core_stops(stops)
        isempty(core) && return Set{Tuple{Int,Int}}()

        idxs = findall(==(y), core)
        isempty(idxs) && return Set{Tuple{Int,Int}}()

        n = size(shortest_old, 1)
        
        pairs_to_check = Set{Tuple{Int,Int}}()
        m = length(core)

        for idx in idxs
            prev_idx = (idx == 1) ? m : idx - 1
            next_idx = (idx == m) ? 1 : idx + 1

            prev_stop = core[prev_idx]
            next_stop = core[next_idx]

            for s in core
                s == y && continue
                push!(pairs_to_check, (y, s))
                push!(pairs_to_check, (s, y))
            end

            left_side  = core[1:prev_idx]
            right_side = core[next_idx:end]

            for u in left_side, v in right_side
                u == v && continue
                push!(pairs_to_check, (u, v))
                push!(pairs_to_check, (v, u))
            end

            push!(pairs_to_check, (prev_stop, y))
            push!(pairs_to_check, (y, prev_stop))
            push!(pairs_to_check, (next_stop, y))
            push!(pairs_to_check, (y, next_stop))
            push!(pairs_to_check, (prev_stop, next_stop))
            push!(pairs_to_check, (next_stop, prev_stop))
        end

        worsen_edges = Set{Tuple{Int,Int}}()
        lineonly_changed_edges = Set{Tuple{Int,Int}}()

        for (u, v) in pairs_to_check
            u == v && continue

            # Check candidate edge u -> v
            old = directlink_old[u, v]
            new = directlink_new[u, v]

            # Edge detected if disappeared or worsened
            if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
                push!(worsen_edges, (u, v))
            end

            # Edge detected if only the assigned line changes
            if has_lineinfo &&
               isfinite(old) && isfinite(new) &&
               abs(old - new) <= tol &&
               lineinfo_old[u, v] != lineinfo_new[u, v]
                push!(lineonly_changed_edges, (u, v))
            end
        end

        # Convert detected direct-link changes into endpoint-based affected OD pairs
        changed_edges = union(worsen_edges, lineonly_changed_edges)
        return direct_link_endpoint_pairs_symmetric(changed_edges)
    end

    # Non-loop case
    remaining = [s for s in stops if s != y]
    isempty(remaining) && return Set{Tuple{Int,Int}}()

    worsen_edges = Set{Tuple{Int,Int}}()
    lineonly_changed_edges = Set{Tuple{Int,Int}}()

    # Check edges incident to removed terminal y
    @inbounds for s in remaining

        # Check incident edge y -> s
        old = directlink_old[y, s]
        new = directlink_new[y, s]

        # Edge detected if removed or worsened
        if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
            push!(worsen_edges, (y, s))
        end

        # Check incident edge s -> y
        old = directlink_old[s, y]
        new = directlink_new[s, y]

        # Edge detected if removed or worsened
        if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
            push!(worsen_edges, (s, y))
        end

        # Edge detected if only the assigned line changes
        if has_lineinfo
            oldw = directlink_old[y, s]
            neww = directlink_new[y, s]
            if isfinite(oldw) && isfinite(neww) &&
               abs(oldw - neww) <= tol &&
               lineinfo_old[y, s] != lineinfo_new[y, s]
                push!(lineonly_changed_edges, (y, s))
            end

            oldw = directlink_old[s, y]
            neww = directlink_new[s, y]
            if isfinite(oldw) && isfinite(neww) &&
               abs(oldw - neww) <= tol &&
               lineinfo_old[s, y] != lineinfo_new[s, y]
                push!(lineonly_changed_edges, (s, y))
            end
        end
    end

    isempty(worsen_edges) && isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    changed_edges = union(worsen_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)

end

# Line removal 
function detect_line_removal(shortest_old::AbstractMatrix{Float64},
                             directlink_old::AbstractMatrix{Float64},
                             directlink_new::AbstractMatrix{Float64},
                             lineplan_old,
                             L_a::Int,
                             tp::Number;
                             tol::Float64 = 1e-9,
                             path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                             lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                             lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    if L_a < 1 || L_a > size(lineplan_old, 1)
        return Set{Tuple{Int,Int}}()
    end

    n = size(shortest_old, 1)
    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    changed_edges = Set{Tuple{Int,Int}}()

    # Check all direct links in the DLN
    @inbounds for u in 1:n
        for v in 1:n
            u == v && continue

            # Check candidate edge u -> v
            old = directlink_old[u, v]
            new = directlink_new[u, v]

            # Edge detected if disappeared or worsened
            if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
                push!(changed_edges, (u, v))
                continue
            end

            # Edge detected if only the assigned line changes
            if has_lineinfo &&
               isfinite(old) && isfinite(new) &&
               abs(old - new) <= tol &&
               lineinfo_old[u, v] != lineinfo_new[u, v]
                push!(changed_edges, (u, v))
            end
        end
    end

    isempty(changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    return direct_link_endpoint_pairs_symmetric(changed_edges)

end

# Insertion 
function detect_insertion_internal(shortest_old::AbstractMatrix{Float64},
                                   directlink_old::AbstractMatrix{Float64},
                                   directlink_new::AbstractMatrix{Float64},
                                   lineplan_old,
                                   line_id::Int,
                                   n_j::Int,
                                   n_k::Int,
                                   x::Int,
                                   tp::Number;
                                   tol::Float64 = 1e-9,
                                   refine::Bool = true,
                                   path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                   lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                   lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    stops = line_stops(lineplan_old, line_id)
    isempty(stops) && return Set{Tuple{Int,Int}}()

    idx_pair = findfirst(i ->
        (stops[i] == n_j && stops[i+1] == n_k) ||
        (stops[i] == n_k && stops[i+1] == n_j),
        1:length(stops)-1
    )
    idx_pair === nothing && return Set{Tuple{Int,Int}}()

    n = size(shortest_old, 1)
    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    pairs_to_check = Set{Tuple{Int,Int}}()
    
    # Broad detection: check all direct links incident to affected nodes
    if !refine
        affected = unique(vcat(stops, [x]))

        # Check outgoing edges from affected nodes
        @inbounds for u in affected
            for v in 1:n
                u == v && continue
                push!(pairs_to_check, (u, v))
            end
        end

        # Check incoming edges to affected nodes
        @inbounds for v in affected
            for u in 1:n
                u == v && continue
                push!(pairs_to_check, (u, v))
            end
        end

    # Refined detection: check edges involving x and links crossing the insertion position
    else
        left  = stops[1:idx_pair]
        right = stops[idx_pair+1:end]

        # Check incident edges s -> x and x -> s
        @inbounds for s in stops
            s == x && continue
            push!(pairs_to_check, (s, x))
            push!(pairs_to_check, (x, s))
        end

        # Check direct links crossing the inserted position
        @inbounds for u in left, v in right
            u == v && continue
            push!(pairs_to_check, (u, v))
            push!(pairs_to_check, (v, u))
        end
    end

    isempty(pairs_to_check) && return Set{Tuple{Int,Int}}()

    worsen_edges           = Set{Tuple{Int,Int}}()
    improve_edges          = Set{Tuple{Int,Int}}()
    lineonly_changed_edges = Set{Tuple{Int,Int}}()

    @inbounds for (u, v) in pairs_to_check
        u == v && continue

        # Check candidate edge u -> v
        old = directlink_old[u, v]
        new = directlink_new[u, v]

        # Edge detected if removed or worsened
        if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
            push!(worsen_edges, (u, v))
        end

        # Edge detected if newly created or improved
        if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
            push!(improve_edges, (u, v))
        end

        # Edge detected if only the assigned line changes
        if has_lineinfo &&
           isfinite(old) && isfinite(new) &&
           abs(old - new) <= tol &&
           lineinfo_old[u, v] != lineinfo_new[u, v]
            push!(lineonly_changed_edges, (u, v))
        end
    end

    isempty(worsen_edges) &&
    isempty(improve_edges) &&
    isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    changed_edges = union(worsen_edges, improve_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)
end

# Removal 
function detect_removal_internal(shortest_old::AbstractMatrix{Float64},
                                 directlink_old::AbstractMatrix{Float64},
                                 directlink_new::AbstractMatrix{Float64},
                                 lineplan_old,
                                 line_id::Int,
                                 y::Int,
                                 tp::Number;
                                 tol::Float64 = 1e-9,
                                 refine::Bool = true,
                                 path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                 lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                 lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    stops = line_stops(lineplan_old, line_id)
    isempty(stops) && return Set{Tuple{Int,Int}}()

    n = size(shortest_old, 1)
    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    pairs_to_check = Set{Tuple{Int,Int}}()

    # Broad detection: check all direct links between stops on the modified line
    if !refine
        @inbounds for u in stops, v in stops
            u == v && continue
            push!(pairs_to_check, (u, v))
        end

    # Refined detection: check edges involving removed node y and links crossing its former position
    else
        # Locate removed node y
        idx = findfirst(==(y), stops)
        idx === nothing && return Set{Tuple{Int,Int}}()

        left  = stops[1:idx-1]
        right = stops[idx+1:end]

        # Check incident edges y -> s and s -> y
        @inbounds for s in stops
            s == y && continue
            push!(pairs_to_check, (y, s))
            push!(pairs_to_check, (s, y))
        end

        # Check direct links crossing the removed position
        @inbounds for u in left, v in right
            u == v && continue
            push!(pairs_to_check, (u, v))
            push!(pairs_to_check, (v, u))
        end
    end

    isempty(pairs_to_check) && return Set{Tuple{Int,Int}}()

    worsen_edges           = Set{Tuple{Int,Int}}()
    improve_edges          = Set{Tuple{Int,Int}}()
    lineonly_changed_edges = Set{Tuple{Int,Int}}()

    @inbounds for (u, v) in pairs_to_check
        u == v && continue

        # Check candidate edge u -> v
        old = directlink_old[u, v]
        new = directlink_new[u, v]

        # Edge detected if removed or worsened
        if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
            push!(worsen_edges, (u, v))
        end

        # Edge detected if newly created or improved
        if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
            push!(improve_edges, (u, v))
        end

        # Edge detected if only the assigned line changes
        if has_lineinfo &&
           isfinite(old) && isfinite(new) &&
           abs(old - new) <= tol &&
           lineinfo_old[u, v] != lineinfo_new[u, v]
            push!(lineonly_changed_edges, (u, v))
        end
    end

    isempty(worsen_edges) &&
    isempty(improve_edges) &&
    isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    changed_edges = union(worsen_edges, improve_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)
end

# Swap within line
function detect_swap_within_line(shortest_old::AbstractMatrix{Float64},
                                 directlink_old::AbstractMatrix{Float64},
                                 directlink_new::AbstractMatrix{Float64},
                                 lineplan_old,
                                 line_id::Int,
                                 n_j::Int,
                                 n_k::Int,
                                 tp::Number;
                                 tol::Float64 = 1e-9,
                                 refine::Bool = true,
                                 path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                 lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                 lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    stops = line_stops(lineplan_old, line_id)

    if !(n_j in stops && n_k in stops) || n_j == n_k
        return Set{Tuple{Int,Int}}()
    end

    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    # Loop-line case
    if is_loop_line(stops)
        core = loop_core_stops(stops)
        isempty(core) && return Set{Tuple{Int,Int}}()

        if !(n_j in core && n_k in core) || n_j == n_k
            return Set{Tuple{Int,Int}}()
        end

        n = size(shortest_old, 1)
        m = length(core)

        idx_j = findfirst(==(n_j), core)
        idx_k = findfirst(==(n_k), core)
        (idx_j === nothing || idx_k === nothing) && return Set{Tuple{Int,Int}}()

        # Helper function: add direct links crossing a circular cut
        function add_circular_cut_crossing_pairs!(pairs::Set{Tuple{Int,Int}},
                                                  core::Vector{Int},
                                                  idx_x::Int,
                                                  n_j::Int,
                                                  n_k::Int)
            mloc = length(core)
            prev_idx = (idx_x == 1) ? mloc : idx_x - 1
            next_idx = (idx_x == mloc) ? 1 : idx_x + 1

            prev_x = core[prev_idx]
            next_x = core[next_idx]

            left_side  = core[1:prev_idx]
            right_side = core[next_idx:end]

            @inbounds for u in left_side, v in right_side
                u == v && continue

                if !((u == n_j && v == n_k) || (u == n_k && v == n_j))
                    push!(pairs, (u, v))
                end
                if !((v == n_j && u == n_k) || (v == n_k && u == n_j))
                    push!(pairs, (v, u))
                end
            end

            push!(pairs, (prev_x, next_x))
            push!(pairs, (next_x, prev_x))
        end

        pairs_to_check = Set{Tuple{Int,Int}}()

        # Broad detection: check all direct links between nodes on the loop
        if !refine
            @inbounds for u in core, v in core
                u == v && continue
                push!(pairs_to_check, (u, v))
            end

        # Refined detection: check edges involving swapped nodes and links crossing the swapped positions
        else

            # Check incident edges involving n_j and n_k
            @inbounds for s in core
                if s != n_j
                    push!(pairs_to_check, (n_j, s))
                    push!(pairs_to_check, (s, n_j))
                end
                if s != n_k
                    push!(pairs_to_check, (n_k, s))
                    push!(pairs_to_check, (s, n_k))
                end
            end

            # Check direct links crossing the swapped positions
            add_circular_cut_crossing_pairs!(pairs_to_check, core, idx_j, n_j, n_k)
            add_circular_cut_crossing_pairs!(pairs_to_check, core, idx_k, n_j, n_k)
        end

        isempty(pairs_to_check) && return Set{Tuple{Int,Int}}()

        worsen_edges           = Set{Tuple{Int,Int}}()
        improve_edges          = Set{Tuple{Int,Int}}()
        lineonly_changed_edges = Set{Tuple{Int,Int}}()

        @inbounds for (u, v) in pairs_to_check
            u == v && continue

            # Check candidate edge u -> v
            old = directlink_old[u, v]
            new = directlink_new[u, v]

            # Edge detected if removed or worsened
            if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
                push!(worsen_edges, (u, v))
            end

            # Edge detected if newly created or improved
            if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
                push!(improve_edges, (u, v))
            end

            # Edge detected if only the assigned line changes
            if has_lineinfo &&
               isfinite(old) && isfinite(new) &&
               abs(old - new) <= tol &&
               lineinfo_old[u, v] != lineinfo_new[u, v]
                push!(lineonly_changed_edges, (u, v))
            end
        end

        isempty(worsen_edges) &&
        isempty(improve_edges) &&
        isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
        changed_edges = union(worsen_edges, improve_edges, lineonly_changed_edges)
        return direct_link_endpoint_pairs_symmetric(changed_edges)
    end

    # Non-loop case
    # Helper function: add direct links crossing a linear cut
    function add_cut_crossing_pairs!(pairs::Set{Tuple{Int,Int}},
                                     stops::Vector{Int},
                                     left_node::Int,
                                     right_node::Int,
                                     n_j::Int,
                                     n_k::Int)

        idx_left  = findfirst(==(left_node), stops)
        idx_right = findfirst(==(right_node), stops)
        (idx_left === nothing || idx_right === nothing) && return

        lo = min(idx_left, idx_right)
        hi = max(idx_left, idx_right)

        left_set  = stops[1:lo]
        right_set = stops[hi:end]

        @inbounds for u in left_set, v in right_set
            u == v && continue

            if !((u == n_j && v == n_k) || (u == n_k && v == n_j))
                push!(pairs, (u, v))
            end
            if !((v == n_j && u == n_k) || (v == n_k && u == n_j))
                push!(pairs, (v, u))
            end
        end
    end

    pairs_to_check = Set{Tuple{Int,Int}}()

    # Broad detection: check all direct links between stops on the line
    if !refine
        @inbounds for u in stops, v in stops
            u == v && continue
            push!(pairs_to_check, (u, v))
        end

    # Refined detection: check edges involving swapped nodes and links crossing their original positions
    else

        # Check incident edges involving n_j and n_k
        @inbounds for s in stops
            if s != n_j
                push!(pairs_to_check, (n_j, s))
                push!(pairs_to_check, (s, n_j))
            end
            if s != n_k
                push!(pairs_to_check, (n_k, s))
                push!(pairs_to_check, (s, n_k))
            end
        end

        idx_j = findfirst(==(n_j), stops)
        idx_k = findfirst(==(n_k), stops)
        (idx_j === nothing || idx_k === nothing) && return Set{Tuple{Int,Int}}()

        # Check direct links crossing the original positions 
        prev_j = (idx_j > 1) ? stops[idx_j - 1] : nothing
        next_j = (idx_j < length(stops)) ? stops[idx_j + 1] : nothing
        if prev_j !== nothing && next_j !== nothing
            add_cut_crossing_pairs!(pairs_to_check, stops, prev_j, next_j, n_j, n_k)
        end

        # Check direct links crossing the original positions
        prev_k = (idx_k > 1) ? stops[idx_k - 1] : nothing
        next_k = (idx_k < length(stops)) ? stops[idx_k + 1] : nothing
        if prev_k !== nothing && next_k !== nothing
            add_cut_crossing_pairs!(pairs_to_check, stops, prev_k, next_k, n_j, n_k)
        end
    end

    isempty(pairs_to_check) && return Set{Tuple{Int,Int}}()

    worsen_edges           = Set{Tuple{Int,Int}}()
    improve_edges          = Set{Tuple{Int,Int}}()
    lineonly_changed_edges = Set{Tuple{Int,Int}}()

    @inbounds for (u, v) in pairs_to_check
        u == v && continue

        # Check candidate edge u -> v
        old = directlink_old[u, v]
        new = directlink_new[u, v]

        # Edge detected if removed or worsened
        if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
            push!(worsen_edges, (u, v))
        end

        # Edge detected if newly created or improved
        if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
            push!(improve_edges, (u, v))
        end

        # Edge detected if only the assigned line changes
        if has_lineinfo &&
           isfinite(old) && isfinite(new) &&
           abs(old - new) <= tol &&
           lineinfo_old[u, v] != lineinfo_new[u, v]
            push!(lineonly_changed_edges, (u, v))
        end
    end

    isempty(worsen_edges) &&
    isempty(improve_edges) &&
    isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    changed_edges = union(worsen_edges, improve_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)
end

# Swap between lines
function detect_swap_between_lines(shortest_old::AbstractMatrix{Float64},
                                   directlink_old::AbstractMatrix{Float64},
                                   directlink_new::AbstractMatrix{Float64},
                                   lineplan_old,
                                   L_a::Int,
                                   L_b::Int,
                                   n_j::Int,
                                   n_k::Int,
                                   tp::Number;
                                   tol::Float64 = 1e-9,
                                   refine::Bool = true,
                                   twoedge::Bool = true,
                                   path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                   lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                   lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    if L_a == L_b
        return Set{Tuple{Int,Int}}()
    end

    stops_a = line_stops(lineplan_old, L_a)
    stops_b = line_stops(lineplan_old, L_b)

    if isempty(stops_a) || isempty(stops_b)
        return Set{Tuple{Int,Int}}()
    end
    if !(n_j in stops_a) || !(n_k in stops_b)
        return Set{Tuple{Int,Int}}()
    end
    if n_j == n_k
        return Set{Tuple{Int,Int}}()
    end

    n = size(shortest_old, 1)
    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    # Helper function: add direct links crossing the removal position
    function add_removal_cut_pairs_allstops!(pairs::Set{Tuple{Int,Int}},
                                             stops::Vector{Int},
                                             idx_x::Int,
                                             n::Int)

        left_set  = (idx_x > 1) ? stops[1:idx_x-1] : Int[]
        right_set = (idx_x < length(stops)) ? stops[idx_x+1:end] : Int[]

        isempty(left_set) && return
        isempty(right_set) && return

        @inbounds for o in 1:n, u in left_set
            o == u && continue
            push!(pairs, (o, u))
            push!(pairs, (u, o))
        end

        @inbounds for d in 1:n, v in right_set
            d == v && continue
            push!(pairs, (d, v))
            push!(pairs, (v, d))
        end

        @inbounds for u in left_set, v in right_set
            u == v && continue
            push!(pairs, (u, v))
            push!(pairs, (v, u))
        end
    end

    pairs_to_check = Set{Tuple{Int,Int}}()

    # Broad detection: check all direct links in the DLN
    if !refine
        @inbounds for u in 1:n, v in 1:n
            u == v && continue
            push!(pairs_to_check, (u, v))
        end

    # Refined detection: check edges involving swapped nodes and links affected by their removal positions
    else

        # Check incident edges involving n_j and n_k
        @inbounds for s in 1:n
            if s != n_j
                push!(pairs_to_check, (n_j, s))
                push!(pairs_to_check, (s, n_j))
            end
            if s != n_k
                push!(pairs_to_check, (n_k, s))
                push!(pairs_to_check, (s, n_k))
            end
        end

        # Check direct links within line L_a
        @inbounds for u in stops_a, v in stops_a
            u == v && continue
            push!(pairs_to_check, (u, v))
        end

        # Check direct links within line L_b
        @inbounds for u in stops_b, v in stops_b
            u == v && continue
            push!(pairs_to_check, (u, v))
        end

        # Check direct links crossing the removal position of n_j
        idx_j = findfirst(==(n_j), stops_a)
        idx_j === nothing && return Set{Tuple{Int,Int}}()
        add_removal_cut_pairs_allstops!(pairs_to_check, stops_a, idx_j, n)

        # Check direct links crossing the removal position of n_k
        idx_k = findfirst(==(n_k), stops_b)
        idx_k === nothing && return Set{Tuple{Int,Int}}()
        add_removal_cut_pairs_allstops!(pairs_to_check, stops_b, idx_k, n)
    end

    isempty(pairs_to_check) && return Set{Tuple{Int,Int}}()

    worsen_edges           = Set{Tuple{Int,Int}}()
    improve_edges          = Set{Tuple{Int,Int}}()
    lineonly_changed_edges = Set{Tuple{Int,Int}}()

    @inbounds for (u, v) in pairs_to_check
        u == v && continue

        # Check candidate edge u -> v
        old = directlink_old[u, v]
        new = directlink_new[u, v]

        # Edge detected if removed or worsened
        if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
            push!(worsen_edges, (u, v))
        end

        # Edge detected if newly created or improved
        if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
            push!(improve_edges, (u, v))
        end

        # Edge detected if only the assigned line changes
        if has_lineinfo &&
           isfinite(old) && isfinite(new) &&
           abs(old - new) <= tol &&
           lineinfo_old[u, v] != lineinfo_new[u, v]
            push!(lineonly_changed_edges, (u, v))
        end
    end

    isempty(worsen_edges) &&
    isempty(improve_edges) &&
    isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    changed_edges = union(worsen_edges, improve_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)
end

# Transfer within line 
function detect_transfer_within_line(shortest_old::AbstractMatrix{Float64},
                                     directlink_old::AbstractMatrix{Float64},
                                     directlink_new::AbstractMatrix{Float64},
                                     lineplan_old,
                                     line_id::Int,
                                     n_j::Int,
                                     n_k::Int,
                                     n_l::Int,
                                     tp::Number;
                                     tol::Float64 = 1e-9,
                                     refine::Bool = true,
                                     path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                     lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                     lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    stops = line_stops(lineplan_old, line_id)

    if !(n_j in stops && n_k in stops && n_l in stops) || n_k == n_l
        return Set{Tuple{Int,Int}}()
    end

    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    # Helper function: add direct links crossing a cut position
    function add_cut_crossing_pairs!(pairs::Set{Tuple{Int,Int}},
                                     stops::Vector{Int},
                                     left_node::Int,
                                     right_node::Int)
        idx_left  = findfirst(==(left_node), stops)
        idx_right = findfirst(==(right_node), stops)
        (idx_left === nothing || idx_right === nothing) && return

        lo = min(idx_left, idx_right)
        hi = max(idx_left, idx_right)

        left_set  = stops[1:lo]
        right_set = stops[hi:end]

        @inbounds for u in left_set, v in right_set
            u == v && continue
            push!(pairs, (u, v))
            push!(pairs, (v, u))
        end
    end

    pairs_to_check = Set{Tuple{Int,Int}}()

    # Broad detection: check all direct links between stops on the line
    if !refine
        @inbounds for u in stops, v in stops
            u == v && continue
            push!(pairs_to_check, (u, v))
        end

    # Refined detection: check edges involving moved node n_j and links affected by the old and new positions
    else

        # Check incident edges involving moved node n_j
        @inbounds for s in stops
            s == n_j && continue
            push!(pairs_to_check, (n_j, s))
            push!(pairs_to_check, (s, n_j))
        end

        # Locate the original position of n_j
        idx_j = findfirst(==(n_j), stops)
        idx_j === nothing && return Set{Tuple{Int,Int}}()

        old_prev = (idx_j > 1) ? stops[idx_j - 1] : nothing
        old_next = (idx_j < length(stops)) ? stops[idx_j + 1] : nothing

        # Check direct links crossing the original position of n_j
        if old_prev !== nothing && old_next !== nothing
            add_cut_crossing_pairs!(pairs_to_check, stops, old_prev, old_next)
        end

        # Check direct links crossing the new insertion position
        add_cut_crossing_pairs!(pairs_to_check, stops, n_k, n_l)
    end

    isempty(pairs_to_check) && return Set{Tuple{Int,Int}}()

    worsen_edges           = Set{Tuple{Int,Int}}()
    improve_edges          = Set{Tuple{Int,Int}}()
    lineonly_changed_edges = Set{Tuple{Int,Int}}()

    @inbounds for (u, v) in pairs_to_check
        u == v && continue

        # Check candidate edge u -> v
        old = directlink_old[u, v]
        new = directlink_new[u, v]

        # Edge detected if removed or worsened
        if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
            push!(worsen_edges, (u, v))
        end

        # Edge detected if newly created or improved
        if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
            push!(improve_edges, (u, v))
        end

        # Edge detected if only the assigned line changes
        if has_lineinfo &&
           isfinite(old) && isfinite(new) &&
           abs(old - new) <= tol &&
           lineinfo_old[u, v] != lineinfo_new[u, v]
            push!(lineonly_changed_edges, (u, v))
        end
    end

    isempty(worsen_edges) &&
    isempty(improve_edges) &&
    isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    changed_edges = union(worsen_edges, improve_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)
end

# Transfer between lines
function detect_transfer_between_lines(shortest_old::AbstractMatrix{Float64},
                                       directlink_old::AbstractMatrix{Float64},
                                       directlink_new::AbstractMatrix{Float64},
                                       lineplan_old,
                                       L_a::Int,
                                       L_b::Int,
                                       n_j::Int,
                                       n_k::Int,
                                       n_l::Int,
                                       tp::Number;
                                       tol::Float64 = 1e-9,
                                       refine::Bool = true,
                                       path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                       lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                       lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    stops_a = line_stops(lineplan_old, L_a)
    stops_b = line_stops(lineplan_old, L_b)

    if isempty(stops_a) || isempty(stops_b)
        return Set{Tuple{Int,Int}}()
    end
    if !(n_j in stops_a) || (n_j in stops_b)
        return Set{Tuple{Int,Int}}()
    end
    if !(n_k in stops_b && n_l in stops_b) || n_k == n_l
        return Set{Tuple{Int,Int}}()
    end

    idx_pair = findfirst(i ->
        (stops_b[i] == n_k && stops_b[i+1] == n_l) ||
        (stops_b[i] == n_l && stops_b[i+1] == n_k),
        1:length(stops_b)-1
    )
    idx_pair === nothing && return Set{Tuple{Int,Int}}()

    n = size(shortest_old, 1)
    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    # Helper function: add direct links crossing an insertion cut
    function add_cut_crossing_pairs_by_index!(pairs::Set{Tuple{Int,Int}},
                                              stops::Vector{Int},
                                              cut_left_end::Int)
        (cut_left_end < 1 || cut_left_end >= length(stops)) && return
        left_set  = stops[1:cut_left_end]
        right_set = stops[cut_left_end+1:end]

        @inbounds for u in left_set, v in right_set
            u == v && continue
            push!(pairs, (u, v))
            push!(pairs, (v, u))
        end
    end

    # Helper function: add direct links crossing a removal cut
    function add_removal_cut_pairs!(pairs::Set{Tuple{Int,Int}},
                                    stops::Vector{Int},
                                    idx_j::Int)
        left_set  = (idx_j > 1) ? stops[1:idx_j-1] : Int[]
        right_set = (idx_j < length(stops)) ? stops[idx_j+1:end] : Int[]

        @inbounds for u in left_set, v in right_set
            u == v && continue
            push!(pairs, (u, v))
            push!(pairs, (v, u))
        end
    end

    pairs_to_check = Set{Tuple{Int,Int}}()

    # Broad detection: check direct links incident to affected nodes
    if !refine
        affected = unique(vcat(stops_a, stops_b, [n_j]))

        @inbounds for u in affected
            for v in 1:n
                u == v && continue
                push!(pairs_to_check, (u, v))
            end
        end

        @inbounds for v in affected
            for u in 1:n
                u == v && continue
                push!(pairs_to_check, (u, v))
            end
        end

    # Refined detection: check edges involving moved node n_j and links affected by the removal and insertion positions
    else
        involved = unique(vcat(stops_a, stops_b))

        # Check incident edges involving moved node n_j
        @inbounds for s in involved
            s == n_j && continue
            push!(pairs_to_check, (n_j, s))
            push!(pairs_to_check, (s, n_j))
        end

        # Check direct links within source line L_a
        @inbounds for u in stops_a, v in stops_a
            u == v && continue
            push!(pairs_to_check, (u, v))
        end

        # Check direct links within destination line L_b
        @inbounds for u in stops_b, v in stops_b
            u == v && continue
            push!(pairs_to_check, (u, v))
        end

        # Check direct links crossing the removal position on L_a
        idx_j = findfirst(==(n_j), stops_a)
        idx_j === nothing && return Set{Tuple{Int,Int}}()
        add_removal_cut_pairs!(pairs_to_check, stops_a, idx_j)

        # Check direct links crossing the insertion position on L_b
        add_cut_crossing_pairs_by_index!(pairs_to_check, stops_b, idx_pair)
    end

    isempty(pairs_to_check) && return Set{Tuple{Int,Int}}()

    worsen_edges           = Set{Tuple{Int,Int}}()
    improve_edges          = Set{Tuple{Int,Int}}()
    lineonly_changed_edges = Set{Tuple{Int,Int}}()

    @inbounds for (u, v) in pairs_to_check
        u == v && continue

        # Check candidate edge u -> v
        old = directlink_old[u, v]
        new = directlink_new[u, v]

        # Edge detected if removed or worsened
        if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
            push!(worsen_edges, (u, v))
        end

        # Edge detected if newly created or improved
        if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
            push!(improve_edges, (u, v))
        end

        # Edge detected if only the assigned line changes
        if has_lineinfo &&
           isfinite(old) && isfinite(new) &&
           abs(old - new) <= tol &&
           lineinfo_old[u, v] != lineinfo_new[u, v]
            push!(lineonly_changed_edges, (u, v))
        end
    end

    isempty(worsen_edges) &&
    isempty(improve_edges) &&
    isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    changed_edges = union(worsen_edges, improve_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)
end

# Substitution 
function detect_substitution(shortest_old::AbstractMatrix{Float64},
                             directlink_old::AbstractMatrix{Float64},
                             directlink_new::AbstractMatrix{Float64},
                             lineplan_old,
                             line_id::Int,
                             n_j::Int,
                             n_k::Int,
                             tp::Number;
                             tol::Float64 = 1e-9,
                             refine::Bool = true,
                             path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                             lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                             lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    stops = line_stops(lineplan_old, line_id)

    if !(n_j in stops) || n_k == n_j
        return Set{Tuple{Int,Int}}()
    end

    idx_j = findfirst(==(n_j), stops)
    idx_j === nothing && return Set{Tuple{Int,Int}}()

    n = size(shortest_old, 1)
    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    # Helper function: add direct links crossing the substitution position
    function add_substitution_cut_pairs!(pairs::Set{Tuple{Int,Int}},
                                         stops::Vector{Int},
                                         idx_j::Int)
        left_set  = (idx_j > 1) ? stops[1:idx_j-1] : Int[]
        right_set = (idx_j < length(stops)) ? stops[idx_j+1:end] : Int[]

        @inbounds for u in left_set, v in right_set
            u == v && continue
            push!(pairs, (u, v))
            push!(pairs, (v, u))
        end
    end

    pairs_to_check = Set{Tuple{Int,Int}}()

    # Broad detection: check all direct links between nodes on the line
    if !refine
        affected = unique(vcat(stops, [n_k]))

        @inbounds for u in affected
            for v in 1:n
                u == v && continue
                push!(pairs_to_check, (u, v))
            end
        end
        @inbounds for v in affected
            for u in 1:n
                u == v && continue
                push!(pairs_to_check, (u, v))
            end
        end
    
    # Refined detection: check edges involving moved node n_j and links affected by the substitution position
    else
        
        # Check incident edges involving removed node n_j
        @inbounds for s in stops
            s == n_j && continue
            push!(pairs_to_check, (n_j, s))
            push!(pairs_to_check, (s, n_j))
        end

        # Check incident edges involving inserted node n_k
        @inbounds for s in stops
            s == n_k && continue
            push!(pairs_to_check, (n_k, s))
            push!(pairs_to_check, (s, n_k))
        end

        # Check direct links crossing the substitution position
        add_substitution_cut_pairs!(pairs_to_check, stops, idx_j)
    end

    isempty(pairs_to_check) && return Set{Tuple{Int,Int}}()

    worsen_edges           = Set{Tuple{Int,Int}}()
    improve_edges          = Set{Tuple{Int,Int}}()
    lineonly_changed_edges = Set{Tuple{Int,Int}}()

    @inbounds for (u, v) in pairs_to_check
        u == v && continue

        # Check candidate edge u -> v
        old = directlink_old[u, v]
        new = directlink_new[u, v]

        # Edge detected if removed or worsened
        if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
            push!(worsen_edges, (u, v))
        end

        # Edge detected if newly created or improved
        if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
            push!(improve_edges, (u, v))
        end

        # Edge detected if only the assigned line changes
        if has_lineinfo &&
           isfinite(old) && isfinite(new) &&
           abs(old - new) <= tol &&
           lineinfo_old[u, v] != lineinfo_new[u, v]
            push!(lineonly_changed_edges, (u, v))
        end
    end

    isempty(worsen_edges) &&
    isempty(improve_edges) &&
    isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    changed_edges = union(worsen_edges, improve_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)
end

# Partial segment swap 
function detect_partial_segment_swap(shortest_old::AbstractMatrix{Float64},
                                     directlink_old::AbstractMatrix{Float64},
                                     directlink_new::AbstractMatrix{Float64},
                                     lineplan_old,
                                     L_a::Int,
                                     n_j::Int,
                                     n_k::Int,
                                     L_b::Int,
                                     n_l::Int,
                                     n_m::Int,
                                     tp::Number;
                                     tol::Float64 = 1e-9,
                                     refine::Bool = true,
                                     path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                     lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                     lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    if L_a == L_b
        return Set{Tuple{Int,Int}}()
    end

    stops_a = line_stops(lineplan_old, L_a)
    stops_b = line_stops(lineplan_old, L_b)

    if isempty(stops_a) || isempty(stops_b)
        return Set{Tuple{Int,Int}}()
    end

    if !(n_j in stops_a && n_k in stops_a) || n_j == n_k
        return Set{Tuple{Int,Int}}()
    end
    if !(n_l in stops_b && n_m in stops_b) || n_l == n_m
        return Set{Tuple{Int,Int}}()
    end

    # Locate swapped segment boundaries on both lines
    idx_j = findfirst(==(n_j), stops_a)
    idx_k = findfirst(==(n_k), stops_a)
    idx_l = findfirst(==(n_l), stops_b)
    idx_m = findfirst(==(n_m), stops_b)

    (idx_j === nothing || idx_k === nothing || idx_l === nothing || idx_m === nothing) &&
        return Set{Tuple{Int,Int}}()

    lo_a, hi_a = min(idx_j, idx_k), max(idx_j, idx_k)
    lo_b, hi_b = min(idx_l, idx_m), max(idx_l, idx_m)

    # Extract swapped segments
    seg_a = stops_a[lo_a:hi_a]
    seg_b = stops_b[lo_b:hi_b]

    seg_a_set = Set(seg_a)
    seg_b_set = Set(seg_b)

    # Remaining nodes outside the swapped segments
    out_a = [s for s in stops_a if !(s in seg_a_set)]
    out_b = [s for s in stops_b if !(s in seg_b_set)]

    n = size(shortest_old, 1)
    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    pairs_to_check = Set{Tuple{Int,Int}}()

    # Broad detection: check all direct links involving any node on the affected lines
    if !refine
        affected = unique(vcat(stops_a, stops_b))

        @inbounds for u in affected
            for v in 1:n
                u == v && continue
                push!(pairs_to_check, (u, v))
            end
        end

        @inbounds for v in affected
            for u in 1:n
                u == v && continue
                push!(pairs_to_check, (u, v))
            end
        end

    # Refined detection: check edges involving swapped nodes and links affected by the removal and insertion positions
    else
        seg_union = unique(vcat(seg_a, seg_b))
        involved  = unique(vcat(stops_a, stops_b))

        # Check incident edges involving nodes in the swapped segments
        @inbounds for u in seg_union
            for v in involved
                u == v && continue
                push!(pairs_to_check, (u, v))
                push!(pairs_to_check, (v, u))
            end
        end

        # Check direct links within line L_a
        @inbounds for u in stops_a, v in stops_a
            u == v && continue
            push!(pairs_to_check, (u, v))
        end

        # Check direct links within line L_b
        @inbounds for u in stops_b, v in stops_b
            u == v && continue
            push!(pairs_to_check, (u, v))
        end

        # Check direct links crossing the swapped segment boundary on L_a
        @inbounds for u in seg_a, v in out_a
            u == v && continue
            push!(pairs_to_check, (u, v))
            push!(pairs_to_check, (v, u))
        end

        # Check direct links crossing the swapped segment boundary on L_b
        @inbounds for u in seg_b, v in out_b
            u == v && continue
            push!(pairs_to_check, (u, v))
            push!(pairs_to_check, (v, u))
        end
    end

    isempty(pairs_to_check) && return Set{Tuple{Int,Int}}()

    worsen_edges           = Set{Tuple{Int,Int}}()
    improve_edges          = Set{Tuple{Int,Int}}()
    lineonly_changed_edges = Set{Tuple{Int,Int}}()

    @inbounds for (u, v) in pairs_to_check
        u == v && continue

        # Check candidate edge u -> v
        old = directlink_old[u, v]
        new = directlink_new[u, v]

        # Edge detected if removed or worsened
        if (isfinite(old) && !isfinite(new)) ||
           (isfinite(old) && isfinite(new) && new > old + tol)
            push!(worsen_edges, (u, v))
        end

        # Edge detected if newly created or improved
        if (!isfinite(old) && isfinite(new)) ||
           (isfinite(old) && isfinite(new) && new < old - tol)
            push!(improve_edges, (u, v))
        end

        # Edge detected if only the assigned line changes
        if has_lineinfo &&
           isfinite(old) && isfinite(new) &&
           abs(old - new) <= tol &&
           lineinfo_old[u, v] != lineinfo_new[u, v]
            push!(lineonly_changed_edges, (u, v))
        end
    end

    isempty(worsen_edges) &&
    isempty(improve_edges) &&
    isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    changed_edges = union(worsen_edges, improve_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)
end

# Reversal 
function detect_reversal(shortest_old::AbstractMatrix{Float64},
                         directlink_old::AbstractMatrix{Float64},
                         directlink_new::AbstractMatrix{Float64},
                         lineplan_old,
                         line_id::Int,
                         tp::Number;
                         tol::Float64 = 1e-9,
                         path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                         lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                         lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    stops = line_stops(lineplan_old, line_id)
    isempty(stops) && return Set{Tuple{Int,Int}}()

    n = size(shortest_old, 1)
    shortest_tp = shortest_to_tp_space(shortest_old, tp)
    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    worsen_edges           = Set{Tuple{Int,Int}}()   # disappeared or increased
    improve_edges          = Set{Tuple{Int,Int}}()   # new or decreased
    lineonly_changed_edges = Set{Tuple{Int,Int}}()   # same cost, different line attribution

    # Check all direct links between stops on the reversed line
    @inbounds for u in stops
        for v in stops
            u == v && continue

            # Check candidate edge u -> v
            old = directlink_old[u, v]
            new = directlink_new[u, v]

            # Edge detected if removed or worsened
            if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
                push!(worsen_edges, (u, v))
            end

            # Edge detected if newly created or improved
            if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
                push!(improve_edges, (u, v))
            end

            # Edge detected if only the assigned line changes
            if has_lineinfo &&
               isfinite(old) && isfinite(new) &&
               abs(old - new) <= tol &&
               lineinfo_old[u, v] != lineinfo_new[u, v]
                push!(lineonly_changed_edges, (u, v))
            end
        end
    end

    isempty(worsen_edges) &&
    isempty(improve_edges) &&
    isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    changed_edges = union(worsen_edges, improve_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)
end

# Segment reversal 
function detect_segment_reversal(shortest_old::AbstractMatrix{Float64},
                                 directlink_old::AbstractMatrix{Float64},
                                 directlink_new::AbstractMatrix{Float64},
                                 lineplan_old,
                                 line_id::Int,
                                 n_i::Int,
                                 n_j::Int,
                                 tp::Number;
                                 tol::Float64 = 1e-9,
                                 refine::Bool = true,
                                 path_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                 lineinfo_old::Union{Nothing,AbstractMatrix{Int}} = nothing,
                                 lineinfo_new::Union{Nothing,AbstractMatrix{Int}} = nothing)

    stops = line_stops(lineplan_old, line_id)

    if !(n_i in stops && n_j in stops) || n_i == n_j
        return Set{Tuple{Int,Int}}()
    end

    # Locate segment boundaries on the line
    idx_i = findfirst(==(n_i), stops)
    idx_j = findfirst(==(n_j), stops)
    (idx_i === nothing || idx_j === nothing) && return Set{Tuple{Int,Int}}()

    has_lineinfo = (lineinfo_old !== nothing && lineinfo_new !== nothing)

    pairs_to_check = Set{Tuple{Int,Int}}()

    # Check all direct links between stops on the affected line
    @inbounds for u in stops, v in stops
        u == v && continue
        push!(pairs_to_check, (u, v))
    end

    isempty(pairs_to_check) && return Set{Tuple{Int,Int}}()

    worsen_edges           = Set{Tuple{Int,Int}}()
    improve_edges          = Set{Tuple{Int,Int}}()
    lineonly_changed_edges = Set{Tuple{Int,Int}}()

    @inbounds for (u, v) in pairs_to_check
        u == v && continue

        # Check candidate edge u -> v
        old = directlink_old[u, v]
        new = directlink_new[u, v]

        # Edge detected if removed or worsened
        if (isfinite(old) && !isfinite(new)) || (isfinite(old) && isfinite(new) && new > old + tol)
            push!(worsen_edges, (u, v))
        end

        # Edge detected if improved or newly created
        if (!isfinite(old) && isfinite(new)) || (isfinite(old) && isfinite(new) && new < old - tol)
            push!(improve_edges, (u, v))
        end

        # Edge detected if only the assigned line changes
        if has_lineinfo &&
           isfinite(old) && isfinite(new) &&
           abs(old - new) <= tol &&
           lineinfo_old[u, v] != lineinfo_new[u, v]
            push!(lineonly_changed_edges, (u, v))
        end
    end

    isempty(worsen_edges) &&
    isempty(improve_edges) &&
    isempty(lineonly_changed_edges) && return Set{Tuple{Int,Int}}()

    # Convert detected direct-link changes into endpoint-based affected OD pairs
    changed_edges = union(worsen_edges, improve_edges, lineonly_changed_edges)
    return direct_link_endpoint_pairs_symmetric(changed_edges)
end