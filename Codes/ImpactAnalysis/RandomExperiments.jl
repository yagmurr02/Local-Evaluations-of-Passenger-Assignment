using Statistics, DataFrames, Tables, Random, XLSX

# Result container
struct MoveSummary
    move::String
    applied_move::String
    total_od_pairs::Int
    affected::Int
    added::Int
    removed::Int
    increased::Int
    decreased::Int
    path_changed_same_cost::Int
    line_changed_same_path::Int
end

# FEASIBILITY CHECKING
did_change(lp, lp2) = size(lp) != size(lp2) || any(lp .!= lp2)

# Check if candidate is feasible, and if not return old lineplan
function enforce_feasible_candidate(lineplan_old, candidate, network, tp, OD;
                                    verbose::Bool=false)
    if !did_change(lineplan_old, candidate)
        return lineplan_old, false
    end

    ok, _, _ = is_feasible_candidate(candidate, network, OD, tp)
    return ok ? (candidate, true) : (lineplan_old, false)
end

# OD-CHANGE ASSESSMENT
# Reconstruct stop sequence for OD pair (i,j) using path matrix from Floyd-Warshall
function reconstruct_fw_nodes_local(i::Int, j::Int, path::AbstractMatrix{Int})
    k = path[i, j]
    if k == 0
        return [i, j]
    else
        left  = reconstruct_fw_nodes_local(i, k, path)
        right = reconstruct_fw_nodes_local(k, j, path)
        return vcat(left[1:end-1], right)
    end
end

# Reconstruct stop sequence and line sequence for OD pair (o,d) using path and lineinfo matrices
function reconstruct_path_and_lines(o::Int,
                                    d::Int,
                                    path::AbstractMatrix{Int},
                                    lineinfo::AbstractMatrix{Int})
    stops = reconstruct_fw_nodes_local(o, d, path)

    lines = Int[]
    for k in 1:length(stops)-1
        u = stops[k]
        v = stops[k+1]
        push!(lines, lineinfo[u, v])
    end

    return stops, lines
end

# Assess OD changes between two lineplans, including path and line changes at same cost
function assess_od_changes_updated(shortest_old::AbstractMatrix{Float64},
                                   shortest_new::AbstractMatrix{Float64},
                                   path_old::AbstractMatrix{Int},
                                   path_new::AbstractMatrix{Int},
                                   lineinfo_old::AbstractMatrix{Int},
                                   lineinfo_new::AbstractMatrix{Int};
                                   tol::Float64 = 1e-9)

    n = size(shortest_old, 1)

    total = n * (n - 1)
    affected = 0
    added = 0
    removed = 0
    increased = 0
    decreased = 0
    path_changed_same_cost = 0
    line_changed_same_path = 0

    for o in 1:n, d in 1:n
        o == d && continue

        sold = shortest_old[o, d]
        snew = shortest_new[o, d]

        if isfinite(sold) != isfinite(snew)
            affected += 1
            if !isfinite(sold) && isfinite(snew)
                added += 1
            elseif isfinite(sold) && !isfinite(snew)
                removed += 1
            end
            continue
        end

        if !isfinite(sold) && !isfinite(snew)
            continue
        end

        if abs(sold - snew) > tol
            affected += 1
            if snew > sold
                increased += 1
            else
                decreased += 1
            end
            continue
        end

        stops_old, lines_old = reconstruct_path_and_lines(o, d, path_old, lineinfo_old)
        stops_new, lines_new = reconstruct_path_and_lines(o, d, path_new, lineinfo_new)

        if stops_old != stops_new
            affected += 1
            path_changed_same_cost += 1
            continue
        end

        if lines_old != lines_new
            affected += 1
            line_changed_same_path += 1
            continue
        end
    end

    return (
        total = total,
        affected = affected,
        added = added,
        removed = removed,
        increased = increased,
        decreased = decreased,
        path_changed_same_cost = path_changed_same_cost,
        line_changed_same_path = line_changed_same_path
    )
end

# Evaluate a single move by comparing shortest paths before and after, and summarizing OD changes
function evaluate_move(move_name, applied_description, lineplan, lineplan_new, network, tp)
    shortest_old, lineinfo_old, _, path_old = shortest_paths(network, lineplan, tp)
    shortest_new, lineinfo_new, _, path_new = shortest_paths(network, lineplan_new, tp)

    stats = assess_od_changes_updated(
        shortest_old,
        shortest_new,
        path_old,
        path_new,
        lineinfo_old,
        lineinfo_new
    )

    return MoveSummary(
        move_name,
        applied_description,
        stats.total,
        stats.affected,
        stats.added,
        stats.removed,
        stats.increased,
        stats.decreased,
        stats.path_changed_same_cost,
        stats.line_changed_same_path
    )
end

# RANDOM MOVE GENERATION AND EXPERIMENTATION FOR IMPACT ANALYSIS
# Insertion at Terminal
function all_valid_insertion_at_terminal(lineplan, network)
    valid_moves = Tuple[]

    num_lines = size(lineplan, 1)
    n_nodes   = size(network, 1)

    connected(a, b) = isfinite(network[a, b]) && network[a, b] > 0

    for L in 1:num_lines
        stops = get_stops(lineplan, L)
        isempty(stops) && continue

        s = first(stops)
        t = last(stops)

        for node in 1:n_nodes
            (node in stops) && continue

            connected(node, s) && push!(valid_moves, (L, node, "start"))
            connected(t, node) && push!(valid_moves, (L, node, "end"))
        end
    end

    return valid_moves
end

function random_insertion_at_terminal(lineplan, network)
    moves = all_valid_insertion_at_terminal(lineplan, network)
    isempty(moves) && error("No valid insertion-at-terminal moves found.")
    L, node, pos = rand(moves)
    desc = "Insertion at Terminal (L=$L, node=$node, pos=$pos)"
    return L, node, pos, desc
end

function run_random_experiments_insertion_at_terminal(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Tuple}()

    attempts = 0
    max_attempts = 20N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        L, node, pos, desc = random_insertion_at_terminal(lineplan, network)

        key = (L, node, pos)
        (key in seen) && continue
        push!(seen, key)

        candidate = insertion_at_terminal(lineplan, L, node; position=pos)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Insertion
function all_valid_insertion(lineplan, network)
    valid_moves = Tuple[]

    num_lines = size(lineplan, 1)
    n_nodes   = size(network, 1)

    connected(a, b) = isfinite(network[a, b]) && network[a, b] > 0

    for L in 1:num_lines
        stops = get_stops(lineplan, L)
        Nst = length(stops)

        for i in 1:(Nst-1)
            n_j = stops[i]
            n_k = stops[i+1]

            for n_l in 1:n_nodes
                (n_l in stops) && continue
                (n_l == n_j || n_l == n_k) && continue

                if connected(n_j, n_l) && connected(n_l, n_k)
                    push!(valid_moves, (L, n_j, n_k, n_l))
                end
            end
        end
    end

    return valid_moves
end

function random_insertion(lineplan, network)
    moves = all_valid_insertion(lineplan, network)
    isempty(moves) && error("No valid insertion moves found.")
    L, n_j, n_k, n_l = rand(moves)
    desc = "Insertion (L=$L, n_j=$n_j, n_k=$n_k, n_l=$n_l)"
    return L, n_j, n_k, n_l, desc
end

function run_random_experiments_insertion(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Tuple}()

    attempts = 0
    max_attempts = 20N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        L, n_j, n_k, n_l, desc = random_insertion(lineplan, network)
        key = (L, n_j, n_k, n_l)
        (key in seen) && continue
        push!(seen, key)

        candidate = insertion(lineplan, L, n_j, n_k, n_l, network)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Removal at Terminal
function all_valid_removal_at_terminal(lineplan)
    valid_moves = Tuple[]
    num_lines = size(lineplan, 1)

    for L in 1:num_lines
        stops = get_stops(lineplan, L)
        if length(stops) > 1
            push!(valid_moves, (L, stops[1], "start"))
            push!(valid_moves, (L, stops[end], "end"))
        end
    end

    return valid_moves
end

function random_removal_at_terminal(lineplan)
    moves = all_valid_removal_at_terminal(lineplan)
    isempty(moves) && error("No valid removal-at-terminal moves found.")
    L, node, pos = rand(moves)
    desc = "Removal at Terminal (L=$L, node=$node, pos=$pos)"
    return L, node, pos, desc
end

function run_random_experiments_removal_at_terminal(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Tuple}()

    attempts = 0
    max_attempts = 20N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        L, node, pos, desc = random_removal_at_terminal(lineplan)
        key = (L, node, pos)
        (key in seen) && continue
        push!(seen, key)

        candidate = removal_at_terminal(lineplan, L; position=pos)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Removal
function all_valid_removal(lineplan)
    valid_moves = Tuple[]
    num_lines = size(lineplan, 1)

    for L in 1:num_lines
        stops = get_stops(lineplan, L)
        length(stops) <= 2 && continue
        for n_j in stops[2:end-1]
            push!(valid_moves, (L, n_j))
        end
    end

    return valid_moves
end

function random_removal(lineplan)
    moves = all_valid_removal(lineplan)
    isempty(moves) && error("No valid removal moves found.")
    L, node = rand(moves)
    desc = "Removal (L=$L, node=$node)"
    return L, node, desc
end

function run_random_experiments_removal(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Tuple}()

    attempts = 0
    max_attempts = 20N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        L, node, desc = random_removal(lineplan)
        key = (L, node)
        (key in seen) && continue
        push!(seen, key)

        candidate = removal(lineplan, L, node)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Swap within Line
function all_valid_swap_within_line(lineplan)
    valid_moves = Tuple[]
    num_lines = size(lineplan, 1)

    for L in 1:num_lines
        stops = get_stops(lineplan, L)
        Nst = length(stops)
        for i in 1:(Nst-1)
            for k in (i+1):Nst
                push!(valid_moves, (L, stops[i], stops[k]))
            end
        end
    end

    return valid_moves
end

function random_swap_within_line(lineplan)
    moves = all_valid_swap_within_line(lineplan)
    isempty(moves) && error("No valid swap-within-line moves found.")
    L, n_j, n_k = rand(moves)
    desc = "Swap within Line (L=$L, n_j=$n_j, n_k=$n_k)"
    return L, n_j, n_k, desc
end

function run_random_experiments_swap_within_line(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Tuple}()

    attempts = 0
    max_attempts = 30N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        L, n_j, n_k, desc = random_swap_within_line(lineplan)

        a = min(n_j, n_k)
        b = max(n_j, n_k)
        key = (L, a, b)
        (key in seen) && continue
        push!(seen, key)

        candidate = swap_within_line(lineplan, L, n_j, n_k, network)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Swap between Lines
function all_valid_swap_between_lines(lineplan)
    valid_moves = Tuple[]
    num_lines = size(lineplan, 1)

    for L_a in 1:num_lines
        stops_a = get_stops(lineplan, L_a)

        for L_b in 1:num_lines
            (L_b == L_a) && continue
            stops_b = get_stops(lineplan, L_b)

            for n_j in stops_a, n_k in stops_b
                (n_j == n_k) && continue
                push!(valid_moves, (L_a, L_b, n_j, n_k))
            end
        end
    end

    return valid_moves
end

function random_swap_between_lines(lineplan)
    moves = all_valid_swap_between_lines(lineplan)
    isempty(moves) && error("No valid swap-between-lines moves found.")
    L_a, L_b, n_j, n_k = rand(moves)
    desc = "Swap Between Lines (L_a=$L_a, L_b=$L_b, n_j=$n_j, n_k=$n_k)"
    return L_a, L_b, n_j, n_k, desc
end

function run_random_experiments_swap_between_lines(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Tuple}()

    attempts = 0
    max_attempts = 40N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        L_a, L_b, n_j, n_k, desc = random_swap_between_lines(lineplan)

        key = (min(L_a, L_b), max(L_a, L_b), min(n_j, n_k), max(n_j, n_k))
        (key in seen) && continue
        push!(seen, key)

        candidate = swap_between_lines(lineplan, L_a, L_b, n_j, n_k, network)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Transfer within Line
function all_valid_transfer_within_line(lineplan, network)
    valid_moves = Tuple[]

    num_lines = size(lineplan, 1)
    n_nodes   = size(network, 1)

    connected(a, b) = isfinite(network[a,b]) && network[a,b] > 0

    function path_exists(a, b)
        visited = Set{Int}()
        queue = [a]
        while !isempty(queue)
            cur = popfirst!(queue)
            cur == b && return true
            push!(visited, cur)
            for x in 1:n_nodes
                if connected(cur, x) && !(x in visited)
                    push!(queue, x)
                end
            end
        end
        return false
    end

    for L in 1:num_lines
        stops = get_stops(lineplan, L)
        Nst = length(stops)

        for n_j in stops
            for i in 1:(Nst-1)
                n_k = stops[i]
                n_l = stops[i+1]
                (n_j == n_k || n_j == n_l) && continue

                if path_exists(n_k, n_j) && path_exists(n_j, n_l)
                    push!(valid_moves, (L, n_j, n_k, n_l))
                end
            end
        end
    end

    return valid_moves
end

function random_transfer_within_line(lineplan, network)
    moves = all_valid_transfer_within_line(lineplan, network)
    isempty(moves) && error("No valid transfer-within-line moves found.")
    L, n_j, n_k, n_l = rand(moves)
    desc = "Transfer Within Line (L=$L, n_j=$n_j, n_k=$n_k, n_l=$n_l)"
    return L, n_j, n_k, n_l, desc
end

function run_random_experiments_transfer_within_line(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Tuple}()

    attempts = 0
    max_attempts = 50N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        L, n_j, n_k, n_l, desc = random_transfer_within_line(lineplan, network)
        key = (L, n_j, n_k, n_l)
        (key in seen) && continue
        push!(seen, key)

        candidate = transfer_within_line(lineplan, L, n_j, n_k, n_l, network)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Transfer between Lines
function all_valid_transfer_between_lines(lineplan, network)
    valid_moves = Tuple[]

    num_lines = size(lineplan, 1)
    n_nodes   = size(network, 1)

    connected(a, b) = isfinite(network[a,b]) && network[a,b] > 0

    function path_exists(a, b)
        visited = Set{Int}()
        queue = [a]
        while !isempty(queue)
            cur = popfirst!(queue)
            cur == b && return true
            push!(visited, cur)
            for x in 1:n_nodes
                if connected(cur, x) && !(x in visited)
                    push!(queue, x)
                end
            end
        end
        return false
    end

    for L_a in 1:num_lines
        stops_a = get_stops(lineplan, L_a)

        for L_b in 1:num_lines
            (L_b == L_a) && continue
            stops_b = get_stops(lineplan, L_b)

            for n_j in stops_a
                (n_j in stops_b) && continue

                for i in 1:(length(stops_b)-1)
                    n_k = stops_b[i]
                    n_l = stops_b[i+1]

                    if path_exists(n_k, n_j) && path_exists(n_j, n_l)
                        push!(valid_moves, (L_a, L_b, n_j, n_k, n_l))
                    end
                end
            end
        end
    end

    return valid_moves
end

function random_transfer_between_lines(lineplan, network)
    moves = all_valid_transfer_between_lines(lineplan, network)
    isempty(moves) && error("No valid transfer-between-lines moves found.")
    L_a, L_b, n_j, n_k, n_l = rand(moves)
    desc = "Transfer Between Lines (L_a=$L_a, L_b=$L_b, n_j=$n_j, n_k=$n_k, n_l=$n_l)"
    return L_a, L_b, n_j, n_k, n_l, desc
end

function run_random_experiments_transfer_between_lines(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Tuple}()

    attempts = 0
    max_attempts = 60N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        L_a, L_b, n_j, n_k, n_l, desc = random_transfer_between_lines(lineplan, network)
        key = (L_a, L_b, n_j, n_k, n_l)
        (key in seen) && continue
        push!(seen, key)

        candidate = transfer_between_lines(lineplan, L_a, L_b, n_j, n_k, n_l, network)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Substitution
function all_valid_substitution(lineplan, network)
    valid_moves = Tuple[]

    num_lines = size(lineplan, 1)
    n_nodes   = size(network, 1)

    connected(a,b) = isfinite(network[a,b]) && network[a,b] > 0

    function path_exists(a,b)
        visited = Set{Int}()
        queue = [a]
        while !isempty(queue)
            cur = popfirst!(queue)
            cur == b && return true
            push!(visited, cur)
            for x in 1:n_nodes
                if connected(cur, x) && !(x in visited)
                    push!(queue, x)
                end
            end
        end
        return false
    end

    for L in 1:num_lines
        stops = get_stops(lineplan, L)

        for idx_j in eachindex(stops)
            n_j = stops[idx_j]

            for n_k in 1:n_nodes
                (n_k == n_j) && continue
                (n_k in stops) && continue

                n_left  = idx_j > 1             ? stops[idx_j - 1] : nothing
                n_right = idx_j < length(stops) ? stops[idx_j + 1] : nothing

                ok_left  = (n_left  === nothing) || path_exists(n_left, n_k)
                ok_right = (n_right === nothing) || path_exists(n_k, n_right)

                if ok_left && ok_right
                    push!(valid_moves, (L, n_j, n_k))
                end
            end
        end
    end

    return valid_moves
end

function random_substitution(lineplan, network)
    moves = all_valid_substitution(lineplan, network)
    isempty(moves) && error("No valid substitution moves found.")
    L, n_j, n_k = rand(moves)
    desc = "Substitution (L=$L, n_j=$n_j, n_k=$n_k)"
    return L, n_j, n_k, desc
end

function run_random_experiments_substitution(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Tuple}()

    attempts = 0
    max_attempts = 40N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        L, n_j, n_k, desc = random_substitution(lineplan, network)
        key = (L, n_j, n_k)
        (key in seen) && continue
        push!(seen, key)

        candidate = substitution(lineplan, L, n_j, n_k, network)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Partial Segment Swap
simple_path_ok(stops::Vector{Int}) = (length(unique(stops)) == length(stops))

function enforce_partial_segment_swap_feasible(lineplan_old, candidate, network, tp, OD)
    did_change(lineplan_old, candidate) || return lineplan_old, false

    ok, _, _ = is_feasible_candidate(candidate, network, OD, tp)
    ok || return lineplan_old, false

    for r in 1:size(candidate, 1)
        st = get_stops(candidate, r)
        simple_path_ok(st) || return lineplan_old, false
    end

    return candidate, true
end

function all_valid_partial_segment_swap(lineplan)
    valid_moves = Tuple[]
    num_lines = size(lineplan, 1)

    for L_a in 1:num_lines
        stops_a = get_stops(lineplan, L_a)
        length(stops_a) < 5 && continue

        for L_b in (L_a+1):num_lines
            stops_b = get_stops(lineplan, L_b)
            length(stops_b) < 5 && continue

            for idx_j in 2:(length(stops_a)-3)
                for idx_k in (idx_j+2):(length(stops_a)-1)
                    n_j = stops_a[idx_j]
                    n_k = stops_a[idx_k]

                    for idx_l in 2:(length(stops_b)-3)
                        for idx_m in (idx_l+2):(length(stops_b)-1)
                            n_l = stops_b[idx_l]
                            n_m = stops_b[idx_m]

                            push!(valid_moves, (L_a, n_j, n_k, L_b, n_l, n_m))
                        end
                    end
                end
            end
        end
    end

    return valid_moves
end

function random_partial_segment_swap(lineplan)
    moves = all_valid_partial_segment_swap(lineplan)
    isempty(moves) && error("No valid partial-segment-swap moves found.")
    L_a, n_j, n_k, L_b, n_l, n_m = rand(moves)
    desc = "Partial Segment Swap (L_a=$L_a, n_j=$n_j, n_k=$n_k, L_b=$L_b, n_l=$n_l, n_m=$n_m)"
    return L_a, n_j, n_k, L_b, n_l, n_m, desc
end

function run_random_experiments_partial_segment_swap(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Tuple}()

    attempts = 0
    max_attempts = 80N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        L_a, n_j, n_k, L_b, n_l, n_m, desc = random_partial_segment_swap(lineplan)
        key = (L_a, n_j, n_k, L_b, n_l, n_m)
        (key in seen) && continue
        push!(seen, key)

        candidate = partial_segment_swap(lineplan, L_a, n_j, n_k, L_b, n_l, n_m, network)
        lineplan_new, ok = enforce_partial_segment_swap_feasible(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Reversal
function all_valid_reversal(lineplan)
    valid = Int[]
    num_lines = size(lineplan, 1)
    for L in 1:num_lines
        stops = get_stops(lineplan, L)
        length(stops) >= 2 && push!(valid, L)
    end
    return valid
end

function random_reversal(lineplan)
    choices = all_valid_reversal(lineplan)
    isempty(choices) && error("No valid reversal moves found.")
    L = rand(choices)
    desc = "Reversal (L=$L)"
    return L, desc
end

function run_random_experiments_reversal(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Int}()

    attempts = 0
    max_attempts = 10N

    while length(results) < N && attempts < max_attempts
        attempts += 1
        L, desc = random_reversal(lineplan)
        (L in seen) && continue
        push!(seen, L)

        candidate = reversal(lineplan, L)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Segment Reversal
function all_valid_segment_reversal(lineplan)
    valid = Tuple[]
    num_lines = size(lineplan, 1)

    for L in 1:num_lines
        stops = get_stops(lineplan, L)
        length(stops) < 3 && continue

        for a in eachindex(stops)
            for b in eachindex(stops)
                n_i = stops[a]
                n_j = stops[b]
                (n_i == n_j) && continue

                _, segment = RemoveSegment_between(stops, n_i, n_j)
                !isempty(segment) && push!(valid, (L, n_i, n_j))
            end
        end
    end

    return valid
end

function random_segment_reversal(lineplan)
    moves = all_valid_segment_reversal(lineplan)
    isempty(moves) && error("No valid segment-reversal moves found.")
    L, n_i, n_j = rand(moves)
    desc = "Segment Reversal (L=$L, n_i=$n_i, n_j=$n_j)"
    return L, n_i, n_j, desc
end

function run_random_experiments_segment_reversal(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Tuple}()

    attempts = 0
    max_attempts = 60N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        L, n_i, n_j, desc = random_segment_reversal(lineplan)
        key = (L, n_i, n_j)
        (key in seen) && continue
        push!(seen, key)

        candidate = segment_reversal(lineplan, L, n_i, n_j)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Line Addition
function random_line_segment(network; min_len=8, max_len=8)
    n_nodes = size(network, 1)
    neighbors(i) = [j for j in 1:n_nodes if isfinite(network[i,j]) && network[i,j] > 0]

    start = rand(1:n_nodes)
    path = [start]

    while length(path) < max_len
        neigh = setdiff(neighbors(path[end]), path)
        isempty(neigh) && break
        push!(path, rand(neigh))
    end

    if length(path) < min_len
        return random_line_segment(network; min_len=min_len, max_len=max_len)
    end

    return path
end

function random_line_addition(lineplan, network)
    n_vec = random_line_segment(network)
    desc = "Line Addition (nodes = $n_vec)"
    return n_vec, desc
end

function run_random_experiments_line_addition(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Tuple}()

    attempts = 0
    max_attempts = 80N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        n_vec, desc = random_line_addition(lineplan, network)
        key = Tuple(n_vec)
        (key in seen) && continue
        push!(seen, key)

        candidate = line_addition(lineplan, n_vec, network)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# Line Removal
function all_valid_line_removal(lineplan)
    num_lines = size(lineplan, 1)
    return [L for L in 1:num_lines if !isempty(get_stops(lineplan, L))]
end

function random_line_removal(lineplan, seen)
    moves = all_valid_line_removal(lineplan)
    available = setdiff(moves, collect(seen))
    isempty(available) && return nothing
    return rand(available)
end

function run_random_experiments_line_removal(move_name, N, lineplan, network, tp, OD)
    results = MoveSummary[]
    seen = Set{Int}()

    attempts = 0
    max_attempts = 20N

    while length(results) < N && attempts < max_attempts
        attempts += 1

        L_a = random_line_removal(lineplan, seen)
        L_a === nothing && break
        push!(seen, L_a)

        candidate = line_removal(lineplan, L_a)
        lineplan_new, ok = enforce_feasible_candidate(lineplan, candidate, network, tp, OD)
        ok || continue

        desc = "Line Removal (L_a=$L_a)"
        push!(results, evaluate_move(move_name, desc, lineplan, lineplan_new, network, tp))
    end

    return results
end

# RESULTS TO DATAFRAME AND SUMMARY STATISTICS
function results_to_dataframe(results)
    return DataFrame(
        move = [r.move for r in results],
        applied_move = [r.applied_move for r in results],
        total_od_pairs = [r.total_od_pairs for r in results],
        affected_od_pairs = [r.affected for r in results],
        added_od_pairs = [r.added for r in results],
        removed_od_pairs = [r.removed for r in results],
        increased_tt = [r.increased for r in results],
        decreased_tt = [r.decreased for r in results],
        path_changed_same_cost = [r.path_changed_same_cost for r in results],
        line_changed_same_path = [r.line_changed_same_path for r in results],
    )
end

function summarize_results(df)
    numeric_cols = [
        :affected_od_pairs,
        :added_od_pairs,
        :removed_od_pairs,
        :increased_tt,
        :decreased_tt,
        :path_changed_same_cost,
        :line_changed_same_path
    ]

    summary_df = DataFrame(
        metric = String[],
        column = Symbol[],
        value = Float64[]
    )

    for col in numeric_cols
        push!(summary_df, ("mean", col, mean(df[!, col])))
        push!(summary_df, ("min", col, minimum(df[!, col])))
        push!(summary_df, ("max", col, maximum(df[!, col])))
        push!(summary_df, ("median", col, median(df[!, col])))
    end

    return summary_df
end

function write_sheet!(filename, sheetname, df; overwrite_file=false)
    df_clean = copy(df)
    for col in names(df_clean)
        if eltype(df_clean[!, col]) == Symbol
            df_clean[!, col] = string.(df_clean[!, col])
        end
    end

    coltable = Tables.columntable(df_clean)

    if overwrite_file || !isfile(filename)
        XLSX.openxlsx(filename, mode="w") do xf
            ws = XLSX.addsheet!(xf, sheetname)
            XLSX.writetable!(ws, coltable)
        end
    else
        XLSX.openxlsx(filename, mode="rw") do xf
            ws = XLSX.addsheet!(xf, sheetname)
            XLSX.writetable!(ws, coltable)
        end
    end
end