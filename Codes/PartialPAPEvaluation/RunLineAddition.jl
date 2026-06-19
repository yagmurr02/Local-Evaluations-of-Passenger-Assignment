# Resolve project root
ROOT = abspath(joinpath(@__DIR__, ".."))

if !isfile(joinpath(ROOT, "evaluator.jl"))
    candidate = joinpath(ROOT, "Codes")
    if isfile(joinpath(candidate, "evaluator.jl"))
        ROOT = abspath(candidate)
    end
end

println("Using ROOT = ", ROOT)
cd(ROOT)

# Includes
include(joinpath(ROOT, "evaluator.jl"))
include(joinpath(ROOT, "detection", "Rules.jl"))
include(joinpath(ROOT, "evaluation", "LocalEvaluation.jl"))

if isfile(joinpath(ROOT, "LocalMoves_BEM.jl"))
    include(joinpath(ROOT, "LocalMoves_BEM.jl"))
elseif isfile(joinpath(ROOT, "BasicElementaryMoves.jl"))
    include(joinpath(ROOT, "BasicElementaryMoves.jl"))
else
    error("Could not find move file.")
end

using DelimitedFiles, DataFrames, XLSX, Statistics, DataStructures, Printf, Random

# Excel writer
function write_df_to_xlsx(path::String, sheetname::String, df::DataFrame; overwrite_file=true)
    mode = overwrite_file ? "w" : "rw"

    XLSX.openxlsx(path, mode=mode) do xf
        if sheetname in XLSX.sheetnames(xf)
            delete!(xf, sheetname)
        end
        ws = XLSX.addsheet!(xf, sheetname)

        for (j, name) in enumerate(names(df))
            ws[1, j] = String(name)
        end

        for i in 1:nrow(df)
            for (j, name) in enumerate(names(df))
                val = df[i, name]
                colname = String(name)

                if startswith(colname, "t_") && val isa AbstractFloat
                    if isnan(val)
                        ws[i+1, j] = "NaN"
                    elseif isinf(val)
                        ws[i+1, j] = val > 0 ? "Inf" : "-Inf"
                    else
                        ws[i+1, j] = @sprintf("%.7f", val)
                    end
                else
                    ws[i+1, j] = val
                end
            end
        end
    end
end

# 5-run average (drop 1st and 5th, average 2-4)
function avg_2_to_4(times::Vector{Float64})
    return mean(times[2:4])
end

# Append average row for selected columns
function append_column_averages(df::DataFrame, cols::Vector{Symbol})
    df_out = deepcopy(df)

    avg_row = Dict{Symbol,Any}()

    for c in names(df_out)
        avg_row[Symbol(c)] = missing
    end

    if :n1 in Symbol.(names(df_out))
        avg_row[:n1] = "AVG"
    elseif :line_id in Symbol.(names(df_out))
        avg_row[:line_id] = "AVG"
    end

    for c in cols
        if c in Symbol.(names(df_out))
            vals = [x for x in df_out[!, c] if !ismissing(x)]
            avg_row[c] = isempty(vals) ? missing : mean(vals)
        end
    end

    push!(df_out, avg_row; cols=:union)
    return df_out
end

# Objective contribution helper
function compute_objective_contribution(affected_set,
                                        shortest_old::AbstractMatrix{Float64},
                                        shortest_new::AbstractMatrix{Float64},
                                        demand::AbstractMatrix{<:Real};
                                        tol::Float64 = 1e-9)

    obj = 0.0

    for (i, j) in affected_set
        i == j && continue

        d = demand[i, j]
        abs(d) <= tol && continue

        old = shortest_old[i, j]
        new = shortest_new[i, j]

        if isfinite(old) && isfinite(new)
            obj += d * (new - old)
        elseif !isfinite(old) && !isfinite(new)
            # no contribution
        else
            return Inf
        end
    end

    return obj
end

# Final objective value helpers
function compute_total_objective(shortest::AbstractMatrix{Float64},
                                 demand::AbstractMatrix{<:Real};
                                 tol::Float64 = 1e-9)

    obj = 0.0
    n = size(shortest, 1)

    for i in 1:n, j in 1:n
        i == j && continue

        d = demand[i, j]
        abs(d) <= tol && continue

        s = shortest[i, j]

        if isfinite(s)
            obj += d * s
        else
            return Inf
        end
    end

    return obj
end

# Helper to build full shortest path matrix from local Dijkstra evaluation
function build_full_matrix_from_local(shortest_old::AbstractMatrix{Float64},
                                      shortest_local::AbstractMatrix{Float64},
                                      affected_set)

    shortest_combined = copy(shortest_old)

    for (i, j) in affected_set
        shortest_combined[i, j] = shortest_local[i, j]
    end

    return shortest_combined
end

# Helper to build full shortest path matrix from local FW subset evaluation
function build_full_matrix_from_subset(shortest_old::AbstractMatrix{Float64},
                                       shortest_sub,
                                       nodes_sub,
                                       affected_set)

    shortest_combined = copy(shortest_old)

    if shortest_sub === nothing || nodes_sub === nothing
        return shortest_combined
    end

    node_to_idx = Dict(node => idx for (idx, node) in enumerate(nodes_sub))

    for (i, j) in affected_set
        if haskey(node_to_idx, i) && haskey(node_to_idx, j)
            ii = node_to_idx[i]
            jj = node_to_idx[j]
            shortest_combined[i, j] = shortest_sub[ii, jj]
        end
    end

    return shortest_combined
end

# Helpers to reconstruct paths and lines for debugging
function expand_fw_stops(path::AbstractMatrix{Int}, o::Int, d::Int)
    k = path[o, d]
    if k == 0
        return [o, d]
    else
        left  = expand_fw_stops(path, o, k)
        right = expand_fw_stops(path, k, d)
        return vcat(left[1:end-1], right)
    end
end

# Convert stop sequence to line sequence
function stops_to_lines(stops::Vector{Int}, lineinfo::AbstractMatrix{Int})
    if length(stops) <= 1
        return Int[]
    end
    lines = Int[]
    for i in 1:length(stops)-1
        push!(lines, lineinfo[stops[i], stops[i+1]])
    end
    return lines
end

# Check if OD pair has changed in cost, path, or lines
function od_changed(od::Tuple{Int,Int},
                    shortest_old::AbstractMatrix{Float64},
                    shortest_new::AbstractMatrix{Float64},
                    path_old::AbstractMatrix{Int},
                    path_new::AbstractMatrix{Int},
                    lineinfo_old::AbstractMatrix{Int},
                    lineinfo_new::AbstractMatrix{Int};
                    tol::Float64 = 1e-9)

    o, d = od

    old_cost = shortest_old[o, d]
    new_cost = shortest_new[o, d]

    cost_changed =
        (isfinite(old_cost) != isfinite(new_cost)) ||
        (isfinite(old_cost) && isfinite(new_cost) && abs(old_cost - new_cost) > tol)

    old_stops = isfinite(old_cost) ? expand_fw_stops(path_old, o, d) : Int[]
    new_stops = isfinite(new_cost) ? expand_fw_stops(path_new, o, d) : Int[]

    path_changed = old_stops != new_stops

    old_lines = isfinite(old_cost) ? stops_to_lines(old_stops, lineinfo_old) : Int[]
    new_lines = isfinite(new_cost) ? stops_to_lines(new_stops, lineinfo_new) : Int[]

    line_changed = old_lines != new_lines

    return cost_changed || path_changed || line_changed
end

# Debug print for OD pair details
function print_od_debug(od::Tuple{Int,Int},
                        shortest_old::AbstractMatrix{Float64},
                        shortest_new::AbstractMatrix{Float64},
                        path_old::AbstractMatrix{Int},
                        path_new::AbstractMatrix{Int},
                        lineinfo_old::AbstractMatrix{Int},
                        lineinfo_new::AbstractMatrix{Int})

    o, d = od

    old_cost = shortest_old[o, d]
    new_cost = shortest_new[o, d]

    old_stops = isfinite(old_cost) ? expand_fw_stops(path_old, o, d) : Int[]
    new_stops = isfinite(new_cost) ? expand_fw_stops(path_new, o, d) : Int[]

    old_lines = isfinite(old_cost) ? stops_to_lines(old_stops, lineinfo_old) : Int[]
    new_lines = isfinite(new_cost) ? stops_to_lines(new_stops, lineinfo_new) : Int[]

    println("  OD=(", o, ",", d, ") | cost_old=", old_cost, " cost_new=", new_cost)
    println("    stops_old = ", old_stops)
    println("    lines_old = ", old_lines)
    println("    stops_new = ", new_stops)
    println("    lines_new = ", new_lines)
end

# Helper to check if two nodes are connected in the network
connected(network, a, b) = isfinite(network[a,b]) && network[a,b] > 0

# Helper to get neighbors of a node in the network
function neighbors_of(network, a::Int)
    n = size(network, 1)
    nbrs = Int[]
    for b in 1:n
        if a != b && connected(network, a, b)
            push!(nbrs, b)
        end
    end
    return nbrs
end

# Helper to check if a sequence of stops forms a simple path (no repeated stops)
simple_path_ok(stops::Vector{Int}) = length(unique(stops)) == length(stops)

# Enumerate feasible Line Addition moves
function enumerate_line_addition_moves(lineplan, network, demand, tp;
                                       max_moves=1000,
                                       max_expansions=200_000,
                                       rng_seed=1,
                                       silent_moves=true)

    rng = MersenneTwister(rng_seed)
    n = size(network,1)

    moves = NamedTuple[]
    expansions = Ref(0)

    function dfs(path::Vector{Int})
        if length(moves) >= max_moves || expansions[] >= max_expansions
            return
        end

        expansions[] += 1

        if length(path) == 25
            n_vec = copy(path)

            lp_new = silent_moves ? redirect_stdout(devnull) do
                line_addition(lineplan, n_vec, network; allow_disconnected=false)
            end : line_addition(lineplan, n_vec, network; allow_disconnected=false)

            lp_new == lineplan && return

            feasible, _, _ = is_feasible_candidate(lp_new, network, demand, tp)
            feasible || return

            for r in 1:size(lp_new,1)
                st = Int.(filter(!=(0), lp_new[r,:]))
                simple_path_ok(st) || return
            end

            push!(moves, (n_vec=n_vec, lineplan_new=lp_new))
            return
        end

        last = path[end]
        nbrs = neighbors_of(network, last)
        shuffle!(rng, nbrs)

        for nb in nbrs
            if nb in path
                continue
            end
            push!(path, nb)
            dfs(path)
            pop!(path)

            if length(moves) >= max_moves || expansions[] >= max_expansions
                return
            end
        end
    end

    start_nodes = shuffle(rng, collect(1:n))

    for s in start_nodes
        length(moves) >= max_moves && break
        expansions[] >= max_expansions && break
        dfs([s])
    end

    println("Enumerated line-addition moves: ", length(moves))
    return moves
end

# Main sweep
function run_line_addition_sweep(; tp=5,
                                 network_file="MandlNetwork.txt",
                                 demand_file="MandlDemand.txt",
                                 out_xlsx="LineAddition_Sweep.xlsx",
                                 max_moves=1000)

    network = readdlm(joinpath(ROOT, network_file))
    demand = Float64.(readdlm(joinpath(ROOT, demand_file)))

    lineplan = [
        1 2 3 6 8 10 11 12 0 0;
        1 2 5 4 6 8 10 11 0 0;
        10 7 15 6 3 2 4 12 0 0;
        9 15 8 10 14 13 11 12 0 0
    ]

    shortest_old, lineinfo_old, directlink_old, path_old = shortest_paths(network, lineplan, tp)

    moves = enumerate_line_addition_moves(lineplan, network, demand, tp;
                                          max_moves=max_moves)

    results = DataFrame(
        n1=Int[], n2=Int[], n3=Int[], n4=Int[],
        n5=Int[], n6=Int[], n7=Int[], n8=Int[],

        n_true=Int[],
        n_rule=Int[],
        n_missed=Int[],
        n_extra=Int[],

        obj_contrib_gt=Float64[],
        obj_contrib_rule=Float64[],
        obj_contrib_gap=Float64[],

        obj_full_eval=Float64[],
        obj_local_dijkstra=Float64[],
        obj_local_dijkstra_pq=Float64[],
        obj_local_fw_subset=Float64[],

        obj_gap_dijkstra=Float64[],
        obj_gap_dijkstra_pq=Float64[],
        obj_gap_fw_subset=Float64[],

        n_subset_nodes=Int[],
        subset_bad=Int[],
        subset_total=Int[],

        t_rule=Float64[],
        t_dln=Float64[],
        t_fw_full=Float64[],
        t_dijkstra=Float64[],
        t_dijkstra_pq=Float64[],
        t_fw_subset=Float64[]
    )

    for m in moves

        shortest_new, directlink_new, lineinfo_new, path_new, _, _ =
            full_evaluation_split(network, m.lineplan_new, tp)

        t_dln_runs = Float64[]
        t_fw_full_runs = Float64[]
        for r in 1:5
            _, _, _, _, tdln_tmp, tfw_tmp =
                full_evaluation_split(network, m.lineplan_new, tp)
            push!(t_dln_runs, tdln_tmp)
            push!(t_fw_full_runs, tfw_tmp)
        end
        t_dln = avg_2_to_4(t_dln_runs)
        t_fw_full = avg_2_to_4(t_fw_full_runs)

        A_true = affected_od_groundtruth(shortest_old, shortest_new,
                                         path_old, path_new,
                                         lineinfo_old, lineinfo_new)

        A_rule = detect_line_addition(shortest_old,
                                      directlink_old,
                                      directlink_new,
                                      m.n_vec,
                                      tp;
                                      path_old=path_old,
                                      lineinfo_old=lineinfo_old,
                                      lineinfo_new=lineinfo_new)

        t_rule_runs = Float64[]
        for r in 1:5
            push!(t_rule_runs, @elapsed begin
                detect_line_addition(shortest_old,
                                     directlink_old,
                                     directlink_new,
                                     m.n_vec,
                                     tp;
                                     path_old=path_old,
                                     lineinfo_old=lineinfo_old,
                                     lineinfo_new=lineinfo_new)
            end)
        end
        t_rule = avg_2_to_4(t_rule_runs)

        n_true, n_rule, n_missed, n_extra =
            detection_quality(A_true, A_rule)

        obj_gt = compute_objective_contribution(A_true,
                                                shortest_old,
                                                shortest_new,
                                                demand)

        obj_rule = compute_objective_contribution(A_rule,
                                                  shortest_old,
                                                  shortest_new,
                                                  demand)

        obj_gap = obj_gt - obj_rule
        obj_full_eval = compute_total_objective(shortest_new, demand)

        shortest_dijkstra_local, _, _, _ =
            local_dijkstra_evaluation(network,
                                      m.lineplan_new,
                                      shortest_old,
                                      tp,
                                      A_rule)

        t_dijkstra_runs = Float64[]
        for r in 1:5
            _, _, _, td_tmp =
                local_dijkstra_evaluation(network,
                                          m.lineplan_new,
                                          shortest_old,
                                          tp,
                                          A_rule)
            push!(t_dijkstra_runs, td_tmp)
        end
        t_dijkstra = avg_2_to_4(t_dijkstra_runs)

        shortest_dijkstra_full =
            build_full_matrix_from_local(shortest_old, shortest_dijkstra_local, A_rule)

        obj_local_dijkstra = compute_total_objective(shortest_dijkstra_full, demand)

        shortest_dijkstra_pq_local, _, _, _ =
            local_dijkstra_pq_evaluation(network,
                                         m.lineplan_new,
                                         shortest_old,
                                         tp,
                                         A_rule)

        t_dijkstra_pq_runs = Float64[]
        for r in 1:5
            _, _, _, tpq_tmp =
                local_dijkstra_pq_evaluation(network,
                                             m.lineplan_new,
                                             shortest_old,
                                             tp,
                                             A_rule)
            push!(t_dijkstra_pq_runs, tpq_tmp)
        end
        t_dijkstra_pq = avg_2_to_4(t_dijkstra_pq_runs)

        shortest_dijkstra_pq_full =
            build_full_matrix_from_local(shortest_old, shortest_dijkstra_pq_local, A_rule)

        obj_local_dijkstra_pq = compute_total_objective(shortest_dijkstra_pq_full, demand)

        shortest_sub, nodes_sub, _, _, _ =
            local_fw_subset_evaluation(network,
                                    lineplan,
                                    m.lineplan_new,
                                    nothing,
                                    tp,
                                    A_rule)

        t_fw_subset_runs = Float64[]
        for r in 1:5
            _, _, _, _, tsub_tmp =
                local_fw_subset_evaluation(network,
                                        lineplan,
                                        m.lineplan_new,
                                        nothing,
                                        tp,
                                        A_rule)
            push!(t_fw_subset_runs, tsub_tmp)
        end
        t_fw_subset = avg_2_to_4(t_fw_subset_runs)

        shortest_fw_subset_full =
            build_full_matrix_from_subset(shortest_old, shortest_sub, nodes_sub, A_rule)

        obj_local_fw_subset = compute_total_objective(shortest_fw_subset_full, demand)

        if shortest_sub === nothing || nodes_sub === nothing
            bad, total = 0, 0
        else
            bad, total =
                subset_quality(shortest_new,
                               shortest_sub,
                               nodes_sub,
                               A_rule)
        end

        obj_gap_dijkstra = obj_local_dijkstra - obj_full_eval
        obj_gap_dijkstra_pq = obj_local_dijkstra_pq - obj_full_eval
        obj_gap_fw_subset = obj_local_fw_subset - obj_full_eval

        n_subset_nodes = nodes_sub === nothing ? 0 : length(nodes_sub)

        nv = m.n_vec

        push!(results, (
            nv[1], nv[2], nv[3], nv[4],
            nv[5], nv[6], nv[7], nv[8],

            n_true, n_rule, n_missed, n_extra,

            obj_gt, obj_rule, obj_gap,

            obj_full_eval, obj_local_dijkstra, obj_local_dijkstra_pq, obj_local_fw_subset,

            obj_gap_dijkstra, obj_gap_dijkstra_pq, obj_gap_fw_subset,

            n_subset_nodes, bad, total,

            t_rule, t_dln, t_fw_full, t_dijkstra, t_dijkstra_pq, t_fw_subset
        ))
    end

    results.t_full_total = results.t_dln .+ results.t_fw_full
    results.t_local_dijkstra_total = results.t_dln .+ results.t_rule .+ results.t_dijkstra
    results.t_local_dijkstra_pq_total = results.t_dln .+ results.t_rule .+ results.t_dijkstra_pq
    results.t_local_fw_subset_total = results.t_dln .+ results.t_rule .+ results.t_fw_subset

    avg_cols = [
        :n_missed,
        :n_extra,
        :obj_contrib_gap,
        :obj_gap_dijkstra,
        :obj_gap_dijkstra_pq,
        :obj_gap_fw_subset,
        :subset_bad,
        :subset_total,
        :t_rule,
        :t_dln,
        :t_fw_full,
        :t_dijkstra,
        :t_dijkstra_pq,
        :t_fw_subset,
        :t_full_total,
        :t_local_dijkstra_total,
        :t_local_dijkstra_pq_total,
        :t_local_fw_subset_total
    ]

    results_out = append_column_averages(results, avg_cols)

    out_path = joinpath(ROOT, out_xlsx)
    write_df_to_xlsx(out_path, "results", results_out; overwrite_file=true)
    println("Wrote: ", out_path)
    println("Moves tested: ", nrow(results))

    return results_out
end

# Run
results = run_line_addition_sweep()