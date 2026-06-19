include("BasicElementaryMoves.jl")

# Infrastructure connectivity check
function is_feasible_lineplan(lineplan, infnetwork)
    violations = Tuple{Int, Int, Int}[]
    nblinesplanned = size(lineplan, 1)

    for line_id in 1:nblinesplanned
        stops = filter(!=(0), lineplan[line_id, :])

        for i in 1:length(stops)-1
            stop_i = stops[i]
            stop_j = stops[i+1]

            if !isfinite(infnetwork[stop_i, stop_j]) || infnetwork[stop_i, stop_j] == Inf
                push!(violations, (line_id, stop_i, stop_j))
            end
        end
    end

    return isempty(violations), violations
end

# Non-zero demand service check
function serves_all_demand(lineplan, network, tp, demand)
    shortest, _, _, _ = shortest_paths(network, lineplan, tp)
    n = size(demand, 1)

    for o in 1:n, d in 1:n
        o == d && continue
        if demand[o, d] > 0 && !isfinite(shortest[o, d])
            return false
        end
    end

    return true
end

# Overall feasibility check
function is_feasible_candidate(lineplan, network, demand, tp)
    ok, violations = is_feasible_lineplan(lineplan, network)
    ok || return false, violations, :infrastructure

    if !serves_all_demand(lineplan, network, tp, demand)
        return false, Tuple{Int,Int,Int}[], :service
    end

    return true, Tuple{Int,Int,Int}[], :ok
end

# Insertion at Terminal
function insertion_at_terminal(lineplan, L_a, n_j; position::String) 
    lineplan_new = copy(lineplan) 
    stops = get_stops(lineplan_new, L_a) 
    N = length(stops) 
    
    if position == "start"
        pos_insert = 1
    elseif position == "end"
        pos_insert = N + 1
    end 
    
    # Elementary Move: Insert(n_l, L_a, pos_insert)
    new_stops = Insert_node(stops, n_j, pos_insert)

    set_stops!(lineplan_new, L_a, new_stops) 

    return lineplan_new 
end

# Insertion
function insertion(lineplan, L_a, n_j, n_k, n_l, network)
    lineplan_new = copy(lineplan)

    stops = get_stops(lineplan_new, L_a)
    N = length(stops)

    idx_pair = findfirst(i ->
        (stops[i] == n_j && stops[i+1] == n_k) ||
        (stops[i] == n_k && stops[i+1] == n_j),
        1:N-1
    )

    if idx_pair === nothing
        return lineplan_new
    end

    n_left  = stops[idx_pair]
    n_right = stops[idx_pair + 1]

    connected(a, b) = isfinite(network[a,b]) && network[a,b] > 0
    if !(connected(n_left, n_l) && connected(n_l, n_right))
        return lineplan_new
    end

    pos_insert = idx_pair + 1

    # Elementary Move: Insert(L_a, n_l, pos_insert)
    new_stops = Insert_node(stops, n_l, pos_insert)

    set_stops!(lineplan_new, L_a, new_stops)

    return lineplan_new
end

# Removal at Terminal
function removal_at_terminal(lineplan, L_a; position::String)
    lineplan_new = copy(lineplan)

    stops = get_stops(lineplan_new, L_a)

    if length(stops) <= 2
        return lineplan_new
    end

    if position == "start"
        n_j = stops[1] # terminal at the beginning
        # Elementary Move: Remove(L_a, n_j)
        new_stops = Remove_node(stops, n_j)

    elseif position == "end"
        n_j = stops[end] # terminal at the end
        # Elementary Move: Remove(L_a, n_j)
        new_stops = Remove_node(stops, n_j)
    end

    set_stops!(lineplan_new, L_a, new_stops)

    return lineplan_new
end

# Removal
function removal(lineplan, L_a, n_j)
    lineplan_new = copy(lineplan)

    stops = get_stops(lineplan_new, L_a)

    if !(n_j in stops)
        return lineplan_new
    end

    if length(stops) <= 2
        return lineplan_new
    end

    # Elementary Move: Remove(L_a, n_j)
    new_stops = Remove_node(stops, n_j)

    set_stops!(lineplan_new, L_a, new_stops)

    return lineplan_new
end

# Swap within Line
function swap_within_line(lineplan, L_a, n_j, n_k, network)
    lineplan_new = copy(lineplan)

    stops = get_stops(lineplan_new, L_a)

    if !(n_j in stops) || !(n_k in stops)
        return lineplan_new
    end

    idx_j = findfirst(==(n_j), stops)
    idx_k = findfirst(==(n_k), stops)

    connected(a, b) = isfinite(network[a,b]) && network[a,b] > 0

    function path_exists(a, b)
        visited = Set{Int}()
        queue = [a]
        while !isempty(queue)
            current = popfirst!(queue)
            current == b && return true
            push!(visited, current)
            for x in 1:size(network,1)
                if connected(current, x) && !(x in visited)
                    push!(queue, x)
                end
            end
        end
        return false
    end

    # Elementary Step: Remove(L_a, n_j)
    stops_1 = Remove_node(stops, n_j)

    # Elementary Step: Insert(L_a, n_j, pos_k)
    pos_j_target = idx_k
    if idx_j < idx_k
        pos_j_target -= 1
    end
    stops_2 = Insert_node(stops_1, n_j, pos_j_target)

    # Elementary Step: Remove(L_a, n_k)
    stops_3 = Remove_node(stops_2, n_k)

    # Elementary Step: Insert(L_a, n_k, pos_j)
    pos_k_target = idx_j
    if idx_k < idx_j
        pos_k_target -= 1
    end
    final_stops = Insert_node(stops_3, n_k, pos_k_target)

    for r in 1:length(final_stops)-1
        if !path_exists(final_stops[r], final_stops[r+1])
            return lineplan_new
        end
    end

    set_stops!(lineplan_new, L_a, final_stops)

    return lineplan_new
end

# Swap between Lines
function swap_between_lines(lineplan, L_a, L_b, n_j, n_k, network)
    lineplan_new = copy(lineplan)

    stops_a = get_stops(lineplan_new, L_a)
    stops_b = get_stops(lineplan_new, L_b)

    if L_a == L_b
        return lineplan_new
    end

    if !(n_j in stops_a)
        return lineplan_new
    end
    if !(n_k in stops_b)
        return lineplan_new
    end

    idx_j = findfirst(==(n_j), stops_a)
    idx_k = findfirst(==(n_k), stops_b)

    connected(a, b) = isfinite(network[a,b]) && network[a,b] > 0

    function path_exists(a, b)
        visited = Set{Int}()
        queue = [a]

        while !isempty(queue)
            current = popfirst!(queue)
            current == b && return true

            push!(visited, current)
            for x in 1:size(network,1)
                if connected(current, x) && !(x in visited)
                    push!(queue, x)
                end
            end
        end

        return false
    end

    # Basic Elementary Step: Remove(L_a, n_j)
    stops_a1 = Remove_node(stops_a, n_j)

    # Basic Elementary Step: Remove(L_b, n_k)
    stops_b1 = Remove_node(stops_b, n_k)

    # Basic Elementary Step: Insert(L_a, n_j, pos_k)
    stops_b2 = Insert_node(stops_b1, n_j, idx_k)

    # Basic Elementary Step: Insert(L_a, n_k, pos_j)
    stops_a2 = Insert_node(stops_a1, n_k, idx_j)

    for r in 1:length(stops_a2)-1
        if !path_exists(stops_a2[r], stops_a2[r+1])
            return lineplan_new
        end
    end

    for r in 1:length(stops_b2)-1
        if !path_exists(stops_b2[r], stops_b2[r+1])
            return lineplan_new
        end
    end

    set_stops!(lineplan_new, L_a, stops_a2)
    set_stops!(lineplan_new, L_b, stops_b2)

    return lineplan_new
end

# Transfer within Line
function transfer_within_line(lineplan, L_a, n_j, n_k, n_l, network)
    lineplan_new = copy(lineplan)

    stops = get_stops(lineplan_new, L_a)

    if !(n_j in stops)
        return lineplan_new
    end
    if !(n_k in stops && n_l in stops)
        return lineplan_new
    end
    if n_j == n_k || n_j == n_l
        return lineplan_new
    end

    idx_pair = findfirst(i ->
        (stops[i] == n_k && stops[i+1] == n_l) ||
        (stops[i] == n_l && stops[i+1] == n_k),
        1:length(stops)-1
    )

    if idx_pair === nothing
        return lineplan_new
    end

    pos_insert = idx_pair + 1

    idx_j = findfirst(==(n_j), stops)

    connected(a,b) = isfinite(network[a,b]) && network[a,b] > 0

    function path_exists(a,b)
        visited = Set{Int}()
        queue = [a]

        while !isempty(queue)
            current = popfirst!(queue)
            current == b && return true
            push!(visited, current)

            for x in 1:size(network,1)
                if connected(current, x) && !(x in visited)
                    push!(queue, x)
                end
            end
        end
        return false
    end

    if !(path_exists(n_k, n_j) && path_exists(n_j, n_l))
        return lineplan_new
    end

    # Basic Elementary Step: Remove(L_a, n_j)
    stops_removed = Remove_node(stops, n_j)

    if idx_j < pos_insert
        pos_insert -= 1
    end

    # Basic Elementary Step: Insert(L_a, n_j, pos_insert)
    new_stops = Insert_node(stops_removed, n_j, pos_insert)

    set_stops!(lineplan_new, L_a, new_stops)

    return lineplan_new
end

# Transfer between Lines
function transfer_between_lines(lineplan, L_a, L_b, n_j, n_k, n_l, network)
    lineplan_new = copy(lineplan)

    stops_a = get_stops(lineplan_new, L_a)     # source line
    stops_b = get_stops(lineplan_new, L_b)     # destination line

    if !(n_j in stops_a)
        return lineplan_new
    end
    if n_j in stops_b
        return lineplan_new
    end
    if !(n_k in stops_b && n_l in stops_b)
        return lineplan_new
    end

    idx_pair = findfirst(i ->
        (stops_b[i] == n_k && stops_b[i+1] == n_l) ||
        (stops_b[i] == n_l && stops_b[i+1] == n_k),
        1:length(stops_b)-1
    )

    if idx_pair === nothing
        return lineplan_new
    end

    pos_insert = idx_pair + 1

    connected(a,b) = isfinite(network[a,b]) && network[a,b] > 0

    function path_exists(a,b)
        visited = Set{Int}()
        queue = [a]

        while !isempty(queue)
            current = popfirst!(queue)
            current == b && return true

            push!(visited, current)
            for x in 1:size(network,1)
                if connected(current, x) && !(x in visited)
                    push!(queue, x)
                end
            end
        end

        return false
    end

    if !(path_exists(n_k, n_j) && path_exists(n_j, n_l))
        return lineplan_new
    end

    # Basic Elementary Step: Remove(L_a, n_j)
    stops_a_new = Remove_node(stops_a, n_j)

    # Basic Elementary Step: Insert(L_b, n_j, pos_insert)
    stops_b_new = Insert_node(stops_b, n_j, pos_insert)

    set_stops!(lineplan_new, L_a, stops_a_new)
    set_stops!(lineplan_new, L_b, stops_b_new)

    return lineplan_new
end

# Substitution
function substitution(lineplan, L_a, n_j, n_k, network)
    lineplan_new = copy(lineplan)

    stops = get_stops(lineplan_new, L_a)

    if !(n_j in stops)
        return lineplan_new
    end
    if n_k in stops
        return lineplan_new
    end

    idx_j = findfirst(==(n_j), stops)

    n_left  = idx_j > 1 ?  stops[idx_j - 1] : nothing
    n_right = idx_j < length(stops) ? stops[idx_j + 1] : nothing

    connected(a,b) = isfinite(network[a,b]) && network[a,b] > 0

    function path_exists(a,b)
        visited = Set{Int}()
        queue = [a]

        while !isempty(queue)
            curr = popfirst!(queue)
            curr == b && return true

            push!(visited, curr)
            for x in 1:size(network,1)
                if connected(curr, x) && !(x in visited)
                    push!(queue, x)
                end
            end
        end

        return false
    end

    ok_left  = (n_left  === nothing) || path_exists(n_left, n_k)
    ok_right = (n_right === nothing) || path_exists(n_k, n_right)

    if !(ok_left && ok_right)
        return lineplan_new
    end

    # Basic Elementary Step: Remove(L_a, n_j)
    stops1 = Remove_node(stops, n_j)

    # Basic Elementary Step: Insert(L_a, n_k, idx_j)
    stops2 = Insert_node(stops1, n_k, idx_j)

    set_stops!(lineplan_new, L_a, stops2)

    return lineplan_new
end

# Partial Segment Swap
function partial_segment_swap(lineplan, L_a, n_j, n_k, L_b, n_l, n_m, network)
    lineplan_new = copy(lineplan)

    stops_a = get_stops(lineplan_new, L_a)
    stops_b = get_stops(lineplan_new, L_b)

    if !(n_j in stops_a && n_k in stops_a)
        return lineplan_new
    end
    if !(n_l in stops_b && n_m in stops_b)
        return lineplan_new
    end

    # Basic Elementary Step: RemoveSegment(L_a, n_j, n_k)
    # Basic Elementary Step: RemoveSegment(L_b, n_l, n_m)
    stops_a_new, seg_a = RemoveSegment_between(stops_a, n_j, n_k)
    stops_b_new, seg_b = RemoveSegment_between(stops_b, n_l, n_m)

    connected(a, b) = isfinite(network[a,b]) && network[a,b] > 0

    function path_exists(a, b)
        visited = Set{Int}()
        queue = [a]

        while !isempty(queue)
            curr = popfirst!(queue)
            curr == b && return true

            push!(visited, curr)
            for x in 1:size(network,1)
                if connected(curr, x) && !(x in visited)
                    push!(queue, x)
                end
            end
        end
        return false
    end

    if !isempty(seg_b)
        if !(path_exists(n_j, seg_b[1]) && path_exists(seg_b[end], n_k))
            return lineplan_new
        end
    end

    if !isempty(seg_a)
        if !(path_exists(n_l, seg_a[1]) && path_exists(seg_a[end], n_m))
            return lineplan_new
        end
    end

    # Basic Elementary Step: InsertSegment(L_a, seg_b, after n_j)
    # Basic Elementary Step: InsertSegment(L_b, seg_a, after n_l)
    stops_a_final = InsertSegment_after(stops_a_new, seg_b, n_j)
    stops_b_final = InsertSegment_after(stops_b_new, seg_a, n_l)

    set_stops!(lineplan_new, L_a, stops_a_final)
    set_stops!(lineplan_new, L_b, stops_b_final)

    return lineplan_new
end

# Reversal
function reversal(lineplan, L_a)
    lineplan_new = copy(lineplan)

    stops_a = get_stops(lineplan_new, L_a)

    if length(stops_a) < 2
        return lineplan_new
    end

    # Basic Elementary Step: Reverse(L_a)
    reversed_stops = reverse(stops_a)

    set_stops!(lineplan_new, L_a, reversed_stops)

    return lineplan_new
end

# Segment Reversal
function segment_reversal(lineplan, L_a, n_i, n_j)
    lineplan_new = copy(lineplan)

    stops_a = get_stops(lineplan_new, L_a)

    if !(n_i in stops_a && n_j in stops_a)
        return lineplan_new
    end

    # Basic Elementary Step: RemoveSegment(L_a, n_i, n_j)
    new_stops_a, segment = RemoveSegment_between(stops_a, n_i, n_j)

    if isempty(segment)
        return lineplan_new
    end

    # Basic Elementary Step: Reverse(segment)
    reversed_segment = Reverse_segment(segment)

    # Basic Elementary Step: InsertSegment(L_a, reversed_segment, after n_i)
    new_stops_a = InsertSegment_after(new_stops_a, reversed_segment, n_i)

    set_stops!(lineplan_new, L_a, new_stops_a)

    return lineplan_new
end

# Line Addition
function line_addition(lineplan, n_vec::Vector{Int}, network; allow_disconnected = false)
    lineplan_new = copy(lineplan)

    if length(n_vec) != 8
        return lineplan_new
    end

    connected(a, b) = isfinite(network[a, b]) && network[a, b] > 0

    function path_exists(a, b)
        visited = Set{Int}()
        queue = [a]
        while !isempty(queue)
            current = popfirst!(queue)
            current == b && return true
            push!(visited, current)
            for n in 1:size(network, 1)
                if connected(current, n) && !(n in visited)
                    push!(queue, n)
                end
            end
        end
        return false
    end

    for i in 1:length(n_vec)-1
        n_i, n_j = n_vec[i], n_vec[i+1]
        if !(connected(n_i, n_j) || (allow_disconnected && path_exists(n_i, n_j)))
            return lineplan_new
        end
    end

    # Basic Elementary Step: CreateSegment([n_1, ..., n_k])
    seg_new = CreateSegment(n_vec)

    L_a = size(lineplan_new, 1) + 1
    new_row = zeros(Int, 1, size(lineplan_new, 2))
    lineplan_new = vcat(lineplan_new, new_row)

    # Basic Elementary Step: InsertSegment (write to new line)
    set_stops!(lineplan_new, L_a, seg_new)

    return lineplan_new
end

# Line Removal
function line_removal(lineplan, L_a)

    if L_a < 1 || L_a > size(lineplan, 1)
        return lineplan
    end

    # Basic Elementary Step: CreateSegment(stops of L_a)
    stops = get_stops(lineplan, L_a)

    if isempty(stops)
        return lineplan
    end

    seg_old = CreateSegment(stops)

    # Basic Elementary Step: RemoveSegment (delete line L_a)
    if L_a == 1
        lineplan_new = lineplan[2:end, :]
    elseif L_a == size(lineplan, 1)
        lineplan_new = lineplan[1:end-1, :]
    else
        lineplan_new = vcat(
            lineplan[1:L_a-1, :],
            lineplan[L_a+1:end, :]
        )
    end

    return lineplan_new
end