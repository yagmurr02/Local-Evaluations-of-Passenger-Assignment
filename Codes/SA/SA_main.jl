include(joinpath(@__DIR__, "solution.jl"))
include(joinpath(@__DIR__, "feasibility.jl"))
include(joinpath(@__DIR__, "initialization.jl"))
include(joinpath(@__DIR__, "sa_moves.jl"))
include(joinpath(@__DIR__, "sa_evaluation.jl"))
include(joinpath(@__DIR__, "sa_neighborhood.jl"))
include(joinpath(@__DIR__, "..", "evaluation", "LocalEvaluation.jl"))

using Random, Statistics, Printf, DelimitedFiles, XLSX, DataFrames, Tables

# Numerical cleanup + formatting helpers
function clean_zero(x::Float64; tol::Float64 = 1e-6)
    return abs(x) <= tol ? 0.0 : x
end

fmt0(x::Real) = @sprintf("%.0f", x)
fmt2(x::Real) = @sprintf("%.2f", x)
fmt3(x::Real) = @sprintf("%.3f", x)
fmt4(x::Real) = @sprintf("%.4f", x)
fmt6(x::Real) = @sprintf("%.6f", x)

fmt0_or_missing(x) = ismissing(x) ? "" : fmt0(x)
fmt2_or_missing(x) = ismissing(x) ? "" : fmt2(x)
fmt3_or_missing(x) = ismissing(x) ? "" : fmt3(x)
fmt4_or_missing(x) = ismissing(x) ? "" : fmt4(x)
fmt6_or_missing(x) = ismissing(x) ? "" : fmt6(x)

# Helper to convert sets of OD pairs to sorted string representation for diagnostics
function pairset_to_string(pairs)
    isempty(pairs) && return ""
    vecpairs = sort!(collect(pairs), by = x -> (x[1], x[2]))
    return join(["($(o),$(d))" for (o, d) in vecpairs], ", ")
end

# Feasibility check helper
function lineplan_is_feasible(lineplan, network, demand, n_nodes, tp; min_nodes::Int, max_nodes::Int)
    return is_feasible(lineplan, network, demand, n_nodes;
                       tp = tp,
                       min_nodes = min_nodes,
                       max_nodes = max_nodes,
                       enforce_connected = false,
                       enforce_all_covered = true)
end

# Hybrid evaluation mode helper
function is_hybrid_eval_mode(eval_mode::Symbol)
    return eval_mode in (
        :local_hybrid_basic,
        :local_hybrid_pq,
        :local_hybrid_basic_periodic,
        :local_hybrid_pq_periodic
    )
end

# Summary statistics helper
function safe_stats(x::Vector{Float64}; tol::Float64 = 1e-6)
    if isempty(x)
        return (
            min = missing,
            max = missing,
            avg = missing,
            median = missing
        )
    end

    return (
        min = clean_zero(minimum(x); tol = tol),
        max = clean_zero(maximum(x); tol = tol),
        avg = clean_zero(mean(x); tol = tol),
        median = clean_zero(Statistics.median(x); tol = tol)
    )
end

# Diagnostic helper to analyze the solution in terms of infinite travel times and served demand
function solution_diagnostics(lineplan, network, demand, tp)
    shortest, _, _, _ = shortest_paths(network, lineplan, tp)
    n_nodes = size(network, 1)

    inf_pairs = Set{Tuple{Int,Int}}()
    inf_pairs_with_demand = Set{Tuple{Int,Int}}()

    total_served_passengers = 0.0

    for o in 1:n_nodes, d in 1:n_nodes
        if o == d
            continue
        end

        if isfinite(shortest[o, d])
            total_served_passengers += demand[o, d]
        else
            push!(inf_pairs, (o, d))
            if demand[o, d] > 0
                push!(inf_pairs_with_demand, (o, d))
            end
        end
    end

    return (
        any_inf_od = !isempty(inf_pairs),
        any_inf_od_with_demand = !isempty(inf_pairs_with_demand),
        n_inf_od_with_demand = length(inf_pairs_with_demand),
        inf_od_with_demand_pairs = pairset_to_string(inf_pairs_with_demand),
        total_served_passengers = total_served_passengers
    )
end

# SA acceptance criterion
function accept_move(delta::Float64, T::Float64, rng::AbstractRNG)
    if delta <= 0
        return true
    end
    return rand(rng) < exp(-delta / T)
end

# SA schedule estimation to set T0 and alpha
function estimate_sa_schedule_paper(initial_sol::TransitSolution,
                                    network,
                                    demand,
                                    tp,
                                    n_nodes::Int,
                                    rng::AbstractRNG;
                                    min_nodes::Int = 2,
                                    max_nodes::Int = typemax(Int),
                                    M::Int = 1000,
                                    N::Int = 1000,
                                    P0::Float64 = 0.999,
                                    PN::Float64 = 0.001)

    deltas = Float64[]

    for _ in 1:M
        cand_lineplan, _ = generate_feasible_neighbor(
            initial_sol, network, demand, tp, n_nodes, rng;
            min_nodes = min_nodes,
            max_nodes = max_nodes,
            max_tries = 2000
        )

        cand_lineplan === nothing && continue

        cand = build_solution(cand_lineplan, network, demand, tp)
        delta = cand.obj - initial_sol.obj
        push!(deltas, abs(delta))
    end

    Δave = isempty(deltas) ? 1.0 : mean(deltas)

    T0 = -Δave / log(P0)
    TN = -Δave / log(PN)
    alpha = exp((log(TN) - log(T0)) / N)

    return T0, TN, alpha, Δave
end

# Helper to get the name of a move type
function move_type_name(move)
    if move isa TerminalInsertion
        return "TerminalInsertion"
    elseif move isa TerminalRemoval
        return "TerminalRemoval"
    elseif move isa Reversal
        return "Reversal"
    elseif move === nothing
        return "None"
    else
        return string(typeof(move))
    end
end

# Main SA function
function simulated_annealing(initial_lineplan,
                             network,
                             demand,
                             tp;
                             eval_mode::Symbol = :full,
                             rng::AbstractRNG = MersenneTwister(1),
                             min_nodes::Int = 2,
                             max_nodes::Int = size(initial_lineplan, 2),
                             time_limit_sec::Float64 = 120.0,
                             T0::Union{Nothing,Float64} = nothing,
                             alpha::Union{Nothing,Float64} = nothing,
                             N_cooling::Int = 1000,
                             L_inner::Int = 100,
                             store_history::Bool = false)

    n_nodes = size(network, 1)

    if !lineplan_is_feasible(initial_lineplan, network, demand, n_nodes, tp;
                             min_nodes = min_nodes,
                             max_nodes = max_nodes)
        error("Initial lineplan is infeasible.")
    end

    current = build_solution(initial_lineplan, network, demand, tp)
    best = deepcopy(current)

    local TN, Δave, alpha_use

    if T0 === nothing || alpha === nothing
        T0_est, TN_est, alpha_est, Δave_est = estimate_sa_schedule_paper(
            current, network, demand, tp, n_nodes, rng;
            min_nodes = min_nodes,
            max_nodes = max_nodes,
            M = 1000,
            N = N_cooling,
            P0 = 0.999,
            PN = 0.001
        )

        T = T0 === nothing ? T0_est : T0
        alpha_use = alpha === nothing ? alpha_est : alpha
        TN = TN_est
        Δave = Δave_est
    else
        T = T0
        alpha_use = alpha
        TN = NaN
        Δave = NaN
    end

    accepted = 0
    improved = 0
    iter = 0

    hybrid_evals = 0
    hybrid_fallbacks = 0
    hybrid_direct_accepts = 0

    history = store_history ? NamedTuple[] : nothing

    start_time = time()
    cooling_step = 0

    while (time() - start_time) < time_limit_sec && cooling_step < N_cooling

        for _ in 1:L_inner
            (time() - start_time) >= time_limit_sec && break

            current_lineplan_before = deepcopy(current.lineplan)

            cand_lineplan, move = generate_feasible_neighbor(
                current, network, demand, tp, n_nodes, rng;
                min_nodes = min_nodes,
                max_nodes = max_nodes
            )

            cand_lineplan === nothing && continue

            if !lineplan_is_feasible(cand_lineplan, network, demand, n_nodes, tp;
                                     min_nodes = min_nodes,
                                     max_nodes = max_nodes)
                continue
            end

            iter += 1
            current_obj_before = current.obj

            cand, used_fallback, detection_time = build_solution_with_meta(
                cand_lineplan,
                current,
                move,
                eval_mode,
                network,
                demand,
                tp;
                iteration = iter,
                full_eval_every = 50
            )

            if is_hybrid_eval_mode(eval_mode)
                hybrid_evals += 1
                if used_fallback
                    hybrid_fallbacks += 1
                else
                    hybrid_direct_accepts += 1
                end
            end

            delta = cand.obj - current_obj_before

            accepted_flag = false
            improved_flag = false

            if accept_move(delta, T, rng)
                current = cand
                accepted += 1
                accepted_flag = true

                if current.obj < best.obj
                    best = deepcopy(current)
                    improved += 1
                    improved_flag = true
                end
            end

            if store_history
                push!(history, (
                    iter = iter,
                    move = deepcopy(move),
                    move_type = move_type_name(move),
                    current_lineplan_before = deepcopy(current_lineplan_before),
                    lineplan = deepcopy(cand_lineplan),
                    current_obj_before = current_obj_before,
                    local_obj = cand.obj,
                    delta = delta,
                    accepted = accepted_flag,
                    improved = improved_flag,
                    temperature = T,
                    used_fallback = used_fallback,
                    detection_time = detection_time
                ))
            end
        end

        T *= alpha_use
        cooling_step += 1
    end

    runtime = time() - start_time

    return (
        eval_mode = eval_mode,
        runtime = runtime,
        iterations = iter,
        best_objective = best.obj,
        final_objective = current.obj,
        best_lineplan = deepcopy(best.lineplan),
        final_lineplan = deepcopy(current.lineplan),
        accepted = accepted,
        improved = improved,
        avg_time_per_iteration = runtime / max(iter, 1),
        iterations_per_second = iter / max(runtime, 1e-12),
        hybrid_evals = hybrid_evals,
        hybrid_fallbacks = hybrid_fallbacks,
        hybrid_direct_accepts = hybrid_direct_accepts,
        hybrid_fallback_ratio = hybrid_evals == 0 ? 0.0 : hybrid_fallbacks / hybrid_evals,
        history = history
    )
end

# Offline accuracy evaluation from stored history
function evaluate_accuracy_from_history(history,
                                        network,
                                        demand,
                                        tp;
                                        tol::Float64 = 1e-6)

    tested = 0
    mismatch_count = 0

    for entry in history
        lineplan = entry.lineplan
        local_obj = entry.local_obj

        full_obj, _, _, _, _ = evaluate_full_solution(network, lineplan, demand, tp)

        abs_err = clean_zero(abs(local_obj - full_obj); tol = tol)

        if abs_err > tol
            mismatch_count += 1
        end

        tested += 1
    end

    if tested == 0
        error("No history entries were available for offline accuracy evaluation.")
    end

    return (
        tested = tested,
        mismatch_count = mismatch_count,
        mismatch_rate_pct = clean_zero(100 * mismatch_count / tested; tol = tol)
    )
end

# Detailed iteration log
function build_accuracy_log_from_history(history,
                                         network,
                                         demand,
                                         tp;
                                         tol::Float64 = 1e-6,
                                         detection_tol::Float64 = 1e-9,
                                         print_missed_extra_debug::Bool = false,
                                         max_rows::Union{Int,Nothing} = nothing)

    total_demand = sum(demand)

    rows = DataFrame(
        Iteration = Union{Missing,Int}[],
        MoveType = String[],

        CurrentTTTBefore = String[],
        LocalTTT = String[],
        FullTTT = String[],
        CurrentAVTBefore = String[],
        LocalAVT = String[],
        FullAVT = String[],
        QualityLoss = String[],
        QualityLossPct = String[],

        Accepted = Union{Missing,Bool}[],
        Improved = Union{Missing,Bool}[],
        Temperature = String[],
        UsedFallback = Union{Missing,Bool}[],
        DetectionTimeSec = String[],

        SubsetClaimedImproving = Union{Missing,Bool}[],
        TrulyImproving = Union{Missing,Bool}[],
        WrongHybridDecision = Union{Missing,Bool}[],

        AnyInfOD = Union{Missing,Bool}[],
        AnyInfODWithDemand = Union{Missing,Bool}[],
        NInfODWithDemand = Union{Missing,Int}[],
        InfODWithDemandPairs = String[],
        TotalServedPassengers = String[],

        NTrue = Union{Missing,Int}[],
        NRule = Union{Missing,Int}[],
        NMissed = Union{Missing,Int}[],
        NExtra = Union{Missing,Int}[],
        MissedPairs = String[],
        ExtraPairs = String[]
    )

    quality_loss_vals = Float64[]
    quality_loss_pct_vals = Float64[]
    detection_time_vals = Float64[]

    limited_history = max_rows !== nothing ? history[1:min(length(history), max_rows)] : history

    for entry in limited_history
        full_obj, _, _, _, _ = evaluate_full_solution(network, entry.lineplan, demand, tp)

        local_obj = clean_zero(entry.local_obj; tol = tol)
        full_obj_clean = clean_zero(full_obj; tol = tol)
        current_before = clean_zero(entry.current_obj_before; tol = tol)

        current_avt = isfinite(current_before) ? clean_zero(current_before / total_demand; tol = tol) : Inf
        local_avt   = isfinite(local_obj)      ? clean_zero(local_obj / total_demand; tol = tol)      : Inf
        full_avt    = isfinite(full_obj_clean) ? clean_zero(full_obj_clean / total_demand; tol = tol) : Inf

        quality_loss = clean_zero(local_obj - full_obj_clean; tol = tol)

        quality_loss_pct =
            abs(full_obj_clean) <= tol ? 0.0 :
            clean_zero(100 * abs(quality_loss) / abs(full_obj_clean); tol = tol)

        push!(quality_loss_vals, quality_loss)
        push!(quality_loss_pct_vals, quality_loss_pct)

        detection_time = clean_zero(entry.detection_time; tol = tol)
        push!(detection_time_vals, detection_time)

        subset_claimed_improving = (local_obj - current_before) < -tol
        truly_improving = (full_obj_clean - current_before) < -tol
        wrong_hybrid_decision = subset_claimed_improving && !truly_improving

        shortest_old, lineinfo_old, directlink_old, path_old =
            shortest_paths(network, entry.current_lineplan_before, tp)

        current_sol_before = TransitSolution(
            deepcopy(entry.current_lineplan_before),
            entry.current_obj_before,
            shortest_old,
            directlink_old,
            path_old,
            lineinfo_old
        )

        shortest_new, lineinfo_new, directlink_new, path_new =
            shortest_paths(network, entry.lineplan, tp)

        diag_new = solution_diagnostics(entry.lineplan, network, demand, tp)

        A_true = affected_od_groundtruth(
            shortest_old, shortest_new,
            path_old, path_new,
            lineinfo_old, lineinfo_new;
            tol = detection_tol
        )

        A_rule = detect_affected_ods(
            current_sol_before,
            entry.lineplan,
            entry.move,
            network,
            tp
        )

        if A_rule === nothing
            A_rule = Set{Tuple{Int,Int}}()
        end

        n_true, n_rule, n_missed, n_extra = detection_quality(A_true, A_rule)

        missed = sort(collect(setdiff(A_true, A_rule)))
        extra  = sort(collect(setdiff(A_rule, A_true)))

        if print_missed_extra_debug && (!isempty(missed) || !isempty(extra))
            println("\n============================================================")
            println("Iteration ", entry.iter, " | Move: ", entry.move_type)
            println("============================================================")

            if !isempty(missed)
                println("Missed ODs list = ", missed)
            end

            if !isempty(extra)
                println("Extra ODs list = ", extra)
            end
        end

        push!(rows, (
            entry.iter,
            entry.move_type,

            fmt2(current_before),
            fmt2(local_obj),
            fmt2(full_obj_clean),
            fmt4(current_avt),
            fmt4(local_avt),
            fmt4(full_avt),
            fmt3(quality_loss),
            fmt3(quality_loss_pct),

            entry.accepted,
            entry.improved,
            fmt4(entry.temperature),
            entry.used_fallback,
            fmt6(detection_time),

            subset_claimed_improving,
            truly_improving,
            wrong_hybrid_decision,

            diag_new.any_inf_od,
            diag_new.any_inf_od_with_demand,
            diag_new.n_inf_od_with_demand,
            diag_new.inf_od_with_demand_pairs,
            fmt2(clean_zero(diag_new.total_served_passengers; tol = tol)),

            n_true,
            n_rule,
            n_missed,
            n_extra,
            pairset_to_string(Set(missed)),
            pairset_to_string(Set(extra))
        ))
    end

    avg_quality_loss =
        isempty(quality_loss_vals) ? missing :
        clean_zero(mean(quality_loss_vals); tol = tol)

    avg_quality_loss_pct =
        isempty(quality_loss_pct_vals) ? missing :
        clean_zero(mean(quality_loss_pct_vals); tol = tol)

    avg_detection_time =
        isempty(detection_time_vals) ? missing :
        clean_zero(mean(detection_time_vals); tol = tol)

    push!(rows, (
        missing,
        "Average",

        "",
        "",
        "",
        "",
        "",
        "",
        fmt3_or_missing(avg_quality_loss),
        fmt3_or_missing(avg_quality_loss_pct),

        missing,
        missing,
        "",
        missing,
        fmt6_or_missing(avg_detection_time),

        missing,
        missing,
        missing,

        missing,
        missing,
        missing,
        "",
        "",

        missing,
        missing,
        missing,
        missing,
        "",
        ""
    ))

    return rows
end

# Run comparison for multiple evaluation modes
function run_multi_method_comparison(initial_lineplan,
                                     network,
                                     demand,
                                     tp;
                                     eval_modes::Vector{Symbol},
                                     seed::Int = 1,
                                     n_runs::Int = 10,
                                     min_nodes::Int = 2,
                                     max_nodes::Int = size(initial_lineplan, 2),
                                     time_limit_sec::Float64 = 120.0,
                                     tol::Float64 = 1e-6,
                                     build_iteration_logs::Bool = false,
                                     detection_tol::Float64 = 1e-9,
                                     print_missed_extra_debug::Bool = false,
                                     max_log_rows::Int = 10000)

    n_nodes = size(network, 1)
    total_demand = sum(demand)

    if !lineplan_is_feasible(initial_lineplan, network, demand, n_nodes, tp;
                             min_nodes = min_nodes,
                             max_nodes = max_nodes)
        error("Initial lineplan passed to run_multi_method_comparison is infeasible.")
    end

    initial_sol = build_solution(initial_lineplan, network, demand, tp)
    initial_diag = solution_diagnostics(initial_lineplan, network, demand, tp)

    T0, TN, alpha_paper, Δave = estimate_sa_schedule_paper(
        initial_sol, network, demand, tp, n_nodes, MersenneTwister(seed);
        min_nodes = min_nodes,
        max_nodes = max_nodes,
        M = 1000,
        N = 1000,
        P0 = 0.999,
        PN = 0.001
    )

    sa_rows = DataFrame(
    Method = String[],

    ReportedBestTTT_Min = String[],
    ReportedBestTTT_Max = String[],
    ReportedBestTTT_Avg = String[],
    ReportedBestTTT_Median = String[],

    FullBestTTT_Min = String[],
    FullBestTTT_Max = String[],
    FullBestTTT_Avg = String[],
    FullBestTTT_Median = String[],

    BestGap_Min = String[],
    BestGap_Max = String[],
    BestGap_Avg = String[],
    BestGap_Median = String[],
    BestGapPct_Min = String[],
    BestGapPct_Max = String[],
    BestGapPct_Avg = String[],
    BestGapPct_Median = String[],

    ReportedFinalTTT_Min = String[],
    ReportedFinalTTT_Max = String[],
    ReportedFinalTTT_Avg = String[],
    ReportedFinalTTT_Median = String[],

    FullFinalTTT_Min = String[],
    FullFinalTTT_Max = String[],
    FullFinalTTT_Avg = String[],
    FullFinalTTT_Median = String[],

    FinalGap_Min = String[],
    FinalGap_Max = String[],
    FinalGap_Avg = String[],
    FinalGap_Median = String[],
    FinalGapPct_Min = String[],
    FinalGapPct_Max = String[],
    FinalGapPct_Avg = String[],
        FinalGapPct_Median = String[],

        ReportedBestAVT_Min = String[],
        ReportedBestAVT_Max = String[],
        ReportedBestAVT_Avg = String[],
        ReportedBestAVT_Median = String[],

        FullBestAVT_Min = String[],
        FullBestAVT_Max = String[],
        FullBestAVT_Avg = String[],
        FullBestAVT_Median = String[],

        ReportedFinalAVT_Min = String[],
        ReportedFinalAVT_Max = String[],
        ReportedFinalAVT_Avg = String[],
        ReportedFinalAVT_Median = String[],

        FullFinalAVT_Min = String[],
        FullFinalAVT_Max = String[],
        FullFinalAVT_Avg = String[],
        FullFinalAVT_Median = String[],

        Runtime = String[],

        Iterations_Min = String[],
        Iterations_Max = String[],
        Iterations_Avg = String[],
        Iterations_Median = String[],

        AcceptedMoves_Min = String[],
        AcceptedMoves_Max = String[],
        AcceptedMoves_Avg = String[],
        AcceptedMoves_Median = String[],

        ImprovingMoves_Min = String[],
        ImprovingMoves_Max = String[],
        ImprovingMoves_Avg = String[],
        ImprovingMoves_Median = String[]
    )

    sa_run_rows = DataFrame(
        Method = String[],
        Run = Int[],

        RuntimeSec = String[],
        Iterations = Int[],

        ReportedBestTTT = String[],
        FullBestTTT = String[],
        BestGap = String[],
        BestGapPct = String[],

        ReportedFinalTTT = String[],
        FullFinalTTT = String[],
        FinalGap = String[],
        FinalGapPct = String[],

        ReportedBestAVT = String[],
        FullBestAVT = String[],
        ReportedFinalAVT = String[],
        FullFinalAVT = String[],

        AcceptedMoves = Int[],
        ImprovingMoves = Int[]
    )

    accuracy_rows = DataFrame(
        Method = String[],
        TestedNeighbors = Int[],
        MismatchCount = Int[],
        MismatchRatePct = String[]
    )

    iteration_logs = Dict{String, DataFrame}()
    best_lineplans = Dict{String, Matrix{Int}}()

    println("==================================================")
    println("MULTI-METHOD COMPARISON")
    @printf("Initial TTT                  = %.6f\n", initial_sol.obj)
    @printf("Initial AVT                  = %.6f\n", initial_sol.obj / total_demand)
    @printf("Initial total served demand  = %.2f\n", initial_diag.total_served_passengers)
    println("Initial AnyInfOD             = ", initial_diag.any_inf_od)
    println("Initial AnyInfODWithDemand   = ", initial_diag.any_inf_od_with_demand)
    println("Initial NInfODWithDemand     = ", initial_diag.n_inf_od_with_demand)
    @printf("Initial temperature          = %.6f\n", T0)
    @printf("Cooling factor alpha         = %.8f\n", alpha_paper)
    @printf("Estimated Δave               = %.6f\n", Δave)
    @printf("Time limit per run           = %.2f s\n", time_limit_sec)
    println("Number of reported runs      = ", n_runs)
    println("Methods                      = ", eval_modes)

    for mode in eval_modes
        println("--------------------------------------------------")
        println("Running SA for method: ", mode)

        println("   Warm-up run (not reported)")

        simulated_annealing(
            initial_lineplan,
            network,
            demand,
            tp;
            eval_mode = mode,
            rng = MersenneTwister(seed),
            min_nodes = min_nodes,
            max_nodes = max_nodes,
            time_limit_sec = time_limit_sec,
            T0 = T0,
            alpha = alpha_paper,
            N_cooling = 1000,
            L_inner = 100,
            store_history = false
        )

        runtimes = Float64[]
        iterations_list = Float64[]

        reported_best_ttt_list = Float64[]
        full_best_ttt_list = Float64[]
        best_gap_list = Float64[]
        best_gap_pct_list = Float64[]

        reported_final_ttt_list = Float64[]
        full_final_ttt_list = Float64[]
        final_gap_list = Float64[]
        final_gap_pct_list = Float64[]

        reported_best_avt_list = Float64[]
        full_best_avt_list = Float64[]
        reported_final_avt_list = Float64[]
        full_final_avt_list = Float64[]

        accepted_list = Float64[]
        improved_list = Float64[]

        best_run_obj = Inf
        best_run_lineplan = deepcopy(initial_lineplan)

        total_tested = 0
        total_mismatch = 0

        for run_id in 1:n_runs
            run_seed = seed + run_id

            println("   Run ", run_id, "/", n_runs, " (seed = ", run_seed, ")")

            sa_result = simulated_annealing(
                initial_lineplan,
                network,
                demand,
                tp;
                eval_mode = mode,
                rng = MersenneTwister(run_seed),
                min_nodes = min_nodes,
                max_nodes = max_nodes,
                time_limit_sec = time_limit_sec,
                T0 = T0,
                alpha = alpha_paper,
                N_cooling = 1000,
                L_inner = 100,
                store_history = build_iteration_logs
            )

            reported_best_obj  = clean_zero(sa_result.best_objective; tol = tol)
            reported_final_obj = clean_zero(sa_result.final_objective; tol = tol)

            best_full_obj, _, _, _, _ = evaluate_full_solution(
                network,
                sa_result.best_lineplan,
                demand,
                tp
            )
            best_full_obj = clean_zero(best_full_obj; tol = tol)

            final_full_obj, _, _, _, _ = evaluate_full_solution(
                network,
                sa_result.final_lineplan,
                demand,
                tp
            )
            final_full_obj = clean_zero(final_full_obj; tol = tol)

            best_gap = clean_zero(reported_best_obj - best_full_obj; tol = tol)
            final_gap = clean_zero(reported_final_obj - final_full_obj; tol = tol)

            best_gap_pct =
                abs(best_full_obj) <= tol ? 0.0 :
                clean_zero(100 * abs(best_gap) / abs(best_full_obj); tol = tol)

            final_gap_pct =
                abs(final_full_obj) <= tol ? 0.0 :
                clean_zero(100 * abs(final_gap) / abs(final_full_obj); tol = tol)

            reported_best_avt  = clean_zero(reported_best_obj / total_demand; tol = tol)
            full_best_avt      = clean_zero(best_full_obj / total_demand; tol = tol)
            reported_final_avt = clean_zero(reported_final_obj / total_demand; tol = tol)
            full_final_avt     = clean_zero(final_full_obj / total_demand; tol = tol)

            push!(runtimes, sa_result.runtime)
            push!(iterations_list, Float64(sa_result.iterations))

            push!(reported_best_ttt_list, reported_best_obj)
            push!(full_best_ttt_list, best_full_obj)
            push!(best_gap_list, best_gap)
            push!(best_gap_pct_list, best_gap_pct)

            push!(reported_final_ttt_list, reported_final_obj)
            push!(full_final_ttt_list, final_full_obj)
            push!(final_gap_list, final_gap)
            push!(final_gap_pct_list, final_gap_pct)

            push!(reported_best_avt_list, reported_best_avt)
            push!(full_best_avt_list, full_best_avt)
            push!(reported_final_avt_list, reported_final_avt)
            push!(full_final_avt_list, full_final_avt)

            push!(accepted_list, Float64(sa_result.accepted))
            push!(improved_list, Float64(sa_result.improved))

            push!(sa_run_rows, (
                String(mode),
                run_id,

                fmt4(clean_zero(sa_result.runtime; tol = tol)),
                sa_result.iterations,

                fmt2(reported_best_obj),
                fmt2(best_full_obj),
                fmt3(best_gap),
                fmt3(best_gap_pct),

                fmt2(reported_final_obj),
                fmt2(final_full_obj),
                fmt3(final_gap),
                fmt3(final_gap_pct),

                fmt4(reported_best_avt),
                fmt4(full_best_avt),
                fmt4(reported_final_avt),
                fmt4(full_final_avt),

                sa_result.accepted,
                sa_result.improved
            ))

            if best_full_obj < best_run_obj
                best_run_obj = best_full_obj
                best_run_lineplan = deepcopy(sa_result.best_lineplan)
            end

            if sa_result.history !== nothing
                println("Checking objective accuracy offline for method: ", mode, " | run: ", run_id)

                acc_result = evaluate_accuracy_from_history(
                    sa_result.history,
                    network,
                    demand,
                    tp;
                    tol = tol
                )

                total_tested += acc_result.tested
                total_mismatch += acc_result.mismatch_count
            end

            if build_iteration_logs && run_id == 1 && sa_result.history !== nothing
                iteration_logs[String(mode)] = build_accuracy_log_from_history(
                    sa_result.history,
                    network,
                    demand,
                    tp;
                    tol = tol,
                    detection_tol = detection_tol,
                    print_missed_extra_debug = print_missed_extra_debug,
                    max_rows = max_log_rows
                )
            end
        end

        best_lineplans[String(mode)] = deepcopy(best_run_lineplan)

        runtime_stats            = safe_stats(runtimes; tol = tol)
        iter_stats               = safe_stats(iterations_list; tol = tol)

        reported_best_ttt_stats  = safe_stats(reported_best_ttt_list; tol = tol)
        full_best_ttt_stats      = safe_stats(full_best_ttt_list; tol = tol)
        best_gap_stats           = safe_stats(best_gap_list; tol = tol)
        best_gap_pct_stats       = safe_stats(best_gap_pct_list; tol = tol)

        reported_final_ttt_stats = safe_stats(reported_final_ttt_list; tol = tol)
        full_final_ttt_stats     = safe_stats(full_final_ttt_list; tol = tol)
        final_gap_stats          = safe_stats(final_gap_list; tol = tol)
        final_gap_pct_stats      = safe_stats(final_gap_pct_list; tol = tol)

        reported_best_avt_stats  = safe_stats(reported_best_avt_list; tol = tol)
        full_best_avt_stats      = safe_stats(full_best_avt_list; tol = tol)
        reported_final_avt_stats = safe_stats(reported_final_avt_list; tol = tol)
        full_final_avt_stats     = safe_stats(full_final_avt_list; tol = tol)

        accepted_stats           = safe_stats(accepted_list; tol = tol)
        improved_stats           = safe_stats(improved_list; tol = tol)

        push!(sa_rows, (
            String(mode),

            fmt2_or_missing(reported_best_ttt_stats.min),
            fmt2_or_missing(reported_best_ttt_stats.max),
            fmt2_or_missing(reported_best_ttt_stats.avg),
            fmt2_or_missing(reported_best_ttt_stats.median),

            fmt2_or_missing(full_best_ttt_stats.min),
            fmt2_or_missing(full_best_ttt_stats.max),
            fmt2_or_missing(full_best_ttt_stats.avg),
            fmt2_or_missing(full_best_ttt_stats.median),

            fmt3_or_missing(best_gap_stats.min),
            fmt3_or_missing(best_gap_stats.max),
            fmt3_or_missing(best_gap_stats.avg),
            fmt3_or_missing(best_gap_stats.median),
            fmt3_or_missing(best_gap_pct_stats.min),
            fmt3_or_missing(best_gap_pct_stats.max),
            fmt3_or_missing(best_gap_pct_stats.avg),
            fmt3_or_missing(best_gap_pct_stats.median),

            fmt2_or_missing(reported_final_ttt_stats.min),
            fmt2_or_missing(reported_final_ttt_stats.max),
            fmt2_or_missing(reported_final_ttt_stats.avg),
            fmt2_or_missing(reported_final_ttt_stats.median),

            fmt2_or_missing(full_final_ttt_stats.min),
            fmt2_or_missing(full_final_ttt_stats.max),
            fmt2_or_missing(full_final_ttt_stats.avg),
            fmt2_or_missing(full_final_ttt_stats.median),

            fmt3_or_missing(final_gap_stats.min),
            fmt3_or_missing(final_gap_stats.max),
            fmt3_or_missing(final_gap_stats.avg),
            fmt3_or_missing(final_gap_stats.median),
            fmt3_or_missing(final_gap_pct_stats.min),
            fmt3_or_missing(final_gap_pct_stats.max),
            fmt3_or_missing(final_gap_pct_stats.avg),
            fmt3_or_missing(final_gap_pct_stats.median),

            fmt4_or_missing(reported_best_avt_stats.min),
            fmt4_or_missing(reported_best_avt_stats.max),
            fmt4_or_missing(reported_best_avt_stats.avg),
            fmt4_or_missing(reported_best_avt_stats.median),

            fmt4_or_missing(full_best_avt_stats.min),
            fmt4_or_missing(full_best_avt_stats.max),
            fmt4_or_missing(full_best_avt_stats.avg),
            fmt4_or_missing(full_best_avt_stats.median),

            fmt4_or_missing(reported_final_avt_stats.min),
            fmt4_or_missing(reported_final_avt_stats.max),
            fmt4_or_missing(reported_final_avt_stats.avg),
            fmt4_or_missing(reported_final_avt_stats.median),

            fmt4_or_missing(full_final_avt_stats.min),
            fmt4_or_missing(full_final_avt_stats.max),
            fmt4_or_missing(full_final_avt_stats.avg),
            fmt4_or_missing(full_final_avt_stats.median),

            fmt4_or_missing(runtime_stats.avg),

            fmt2_or_missing(iter_stats.min),
            fmt2_or_missing(iter_stats.max),
            fmt2_or_missing(iter_stats.avg),
            fmt2_or_missing(iter_stats.median),

            fmt2_or_missing(accepted_stats.min),
            fmt2_or_missing(accepted_stats.max),
            fmt2_or_missing(accepted_stats.avg),
            fmt2_or_missing(accepted_stats.median),

            fmt2_or_missing(improved_stats.min),
            fmt2_or_missing(improved_stats.max),
            fmt2_or_missing(improved_stats.avg),
            fmt2_or_missing(improved_stats.median)
        ))

        mismatch_rate_pct =
            total_tested == 0 ? 0.0 : clean_zero(100 * total_mismatch / total_tested; tol = tol)

        push!(accuracy_rows, (
            String(mode),
            total_tested,
            total_mismatch,
            fmt2(mismatch_rate_pct)
        ))
    end

    return (
        initial_ttt = fmt2(clean_zero(initial_sol.obj; tol = tol)),
        initial_avt = fmt4(clean_zero(initial_sol.obj / total_demand; tol = tol)),
        initial_any_inf_od = initial_diag.any_inf_od,
        initial_any_inf_od_with_demand = initial_diag.any_inf_od_with_demand,
        initial_n_inf_od_with_demand = initial_diag.n_inf_od_with_demand,
        initial_total_served_passengers = fmt2(clean_zero(initial_diag.total_served_passengers; tol = tol)),
        initial_temperature = fmt4(clean_zero(T0; tol = tol)),
        initialization_time = nothing,
        initial_lineplan = deepcopy(initial_lineplan),
        best_lineplans = best_lineplans,
        sa_results = sa_rows,
        sa_run_results = sa_run_rows,
        accuracy_results = accuracy_rows,
        iteration_logs = iteration_logs
    )
end

# Export to Excel
function export_comparison_to_excel(filename, results)

    function make_unique_sheet_name(base_name::String, used_names::Set{String})
        max_len = 31
        name = length(base_name) > max_len ? base_name[1:max_len] : base_name

        if !(name in used_names)
            push!(used_names, name)
            return name
        end

        counter = 1
        while true
            suffix = "_" * string(counter)
            trim_len = max_len - length(suffix)
            candidate = (length(base_name) > trim_len ? base_name[1:trim_len] : base_name) * suffix

            if !(candidate in used_names)
                push!(used_names, candidate)
                return candidate
            end

            counter += 1
        end
    end

    XLSX.openxlsx(filename, mode = "w") do xf
        used_sheet_names = Set{String}()

        summary_df = DataFrame(
            Metric = [
                "Initial TTT",
                "Initial AVT",
                "Initial AnyInfOD",
                "Initial AnyInfODWithDemand",
                "Initial NInfODWithDemand",
                "Initial TotalServedPassengers",
                "Initial temperature",
                "Initialization time (sec)"
            ],
            Value = [
                results.initial_ttt,
                results.initial_avt,
                string(results.initial_any_inf_od),
                string(results.initial_any_inf_od_with_demand),
                string(results.initial_n_inf_od_with_demand),
                results.initial_total_served_passengers,
                results.initial_temperature,
                results.initialization_time
            ]
        )

        sheet1_name = make_unique_sheet_name("Summary", used_sheet_names)
        sheet1 = XLSX.addsheet!(xf, sheet1_name)
        XLSX.writetable!(sheet1, Tables.columntable(summary_df); write_columnnames = true)

        sheet2_name = make_unique_sheet_name("SA_Results", used_sheet_names)
        sheet2 = XLSX.addsheet!(xf, sheet2_name)
        XLSX.writetable!(sheet2, Tables.columntable(results.sa_results); write_columnnames = true)

        sheet_runs_name = make_unique_sheet_name("SA_Run_Results", used_sheet_names)
        sheet_runs = XLSX.addsheet!(xf, sheet_runs_name)
        XLSX.writetable!(sheet_runs, Tables.columntable(results.sa_run_results); write_columnnames = true)

        sheet3_name = make_unique_sheet_name("Accuracy_vs_Full", used_sheet_names)
        sheet3 = XLSX.addsheet!(xf, sheet3_name)
        XLSX.writetable!(sheet3, Tables.columntable(results.accuracy_results); write_columnnames = true)

        sheet_init_name = make_unique_sheet_name("Initial_Lineplan", used_sheet_names)
        sheet_init = XLSX.addsheet!(xf, sheet_init_name)
        XLSX.writetable!(
            sheet_init,
            Tables.columntable(DataFrame(results.initial_lineplan, :auto));
            write_columnnames = false
        )

        for (method, lp) in results.best_lineplans
            base_name = "Best_" * replace(method, ":" => "")
            sheet_name = make_unique_sheet_name(base_name, used_sheet_names)

            sheet = XLSX.addsheet!(xf, sheet_name)
            XLSX.writetable!(
                sheet,
                Tables.columntable(DataFrame(lp, :auto));
                write_columnnames = false
            )
        end

        for (method, df) in results.iteration_logs
            base_name = replace(method, ":" => "")
            sheet_name = make_unique_sheet_name(base_name, used_sheet_names)

            sheet = XLSX.addsheet!(xf, sheet_name)
            XLSX.writetable!(sheet, Tables.columntable(df); write_columnnames = true)
        end
    end

    println("==================================================")
    println("Results exported to: ", filename)
end

# Data loading
BASE = normpath(joinpath(@__DIR__, ".."))

network_path  = joinpath(BASE, "TestNetwork.txt")
demand_path   = joinpath(BASE, "TestDemand.txt")
lineplan_path = joinpath(BASE, "TestLineplan.txt")

network         = Float64.(readdlm(network_path))
demand          = Float64.(readdlm(demand_path))
loaded_lineplan = readdlm(lineplan_path, Int)
tp = 5.0

# Build feasible initial lineplan
seed = 1
rng_init = MersenneTwister(seed)

n_nodes = size(network, 1)
n_lines = size(loaded_lineplan, 1)
max_nodes = size(loaded_lineplan, 2)
min_nodes = 2

println("==================================================")
println("INITIALIZATION CHECK")

loaded_feasible = lineplan_is_feasible(loaded_lineplan, network, demand, n_nodes, tp;
                                       min_nodes = min_nodes,
                                       max_nodes = max_nodes)

println("Loaded lineplan feasible? ", loaded_feasible)
println("Constructing a new feasible initial lineplan...")

init_start_time = time()

initial_lineplan = initialize_feasible_lineplan(
    network,
    demand,
    tp,
    n_lines,
    rng_init;
    min_nodes = min_nodes,
    max_nodes = max_nodes,
    outer_attempts = 300,
    max_route_attempts = 300,
    max_cover_steps = 200,
    max_connect_steps = 200,
    verbose = true
)

init_runtime = time() - init_start_time

initialized_feasible = lineplan_is_feasible(initial_lineplan, network, demand, n_nodes, tp;
                                            min_nodes = min_nodes,
                                            max_nodes = max_nodes)

println("Constructed initial lineplan feasible? ", initialized_feasible)

if !initialized_feasible
    error("Constructed initial lineplan is infeasible.")
end

@printf("Initialization time = %.4f seconds\n", init_runtime)

# Methods to compare
eval_modes = [
    :full,
    :local_dijkstra_basic,
    :local_dijkstra_pq,
    :local_floyd_subset,
    :local_floyd_subset_path,
    :local_fw_subset_path_basic,
    :local_fw_subset_path_pq,
    :local_hybrid_basic,
    :local_hybrid_pq,
    :local_threshold_basic,
    :local_threshold_pq,
    :local_floyd_subset_periodic,
    :local_hybrid_basic_periodic,
    :local_hybrid_pq_periodic
]

# Run experiment
results = run_multi_method_comparison(
    initial_lineplan,
    network,
    demand,
    tp;
    eval_modes = eval_modes,
    seed = seed,
    n_runs = 5,
    min_nodes = min_nodes,
    max_nodes = size(initial_lineplan, 2),
    time_limit_sec = 300.0,
    tol = 1e-6,
    build_iteration_logs = false,
    detection_tol = 1e-9,
    print_missed_extra_debug = false,
    max_log_rows = 1000
)

results = merge(results, (initialization_time = fmt4(init_runtime),))

export_comparison_to_excel(
    joinpath(BASE, "evaluation_method_comparison.xlsx"),
    results
)

nothing