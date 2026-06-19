include(joinpath(@__DIR__, "solution.jl"))
include(joinpath(@__DIR__, "feasibility.jl"))
include(joinpath(@__DIR__, "initialization.jl"))

include(joinpath(@__DIR__, "basic_vns_moves.jl"))
include(joinpath(@__DIR__, "basic_vns_neighborhood.jl"))
include(joinpath(@__DIR__, "basic_vns_evaluation.jl"))

include(joinpath(@__DIR__, "..", "evaluation", "LocalEvaluation.jl"))

using Random, Statistics, Printf, DelimitedFiles, XLSX, DataFrames, Tables

# Numerical cleanup + formatting helpers
clean_zero(x::Float64; tol::Float64 = 1e-6) = abs(x) <= tol ? 0.0 : x

fmt2(x::Real) = @sprintf("%.2f", x)
fmt3(x::Real) = @sprintf("%.3f", x)
fmt4(x::Real) = @sprintf("%.4f", x)
fmt6(x::Real) = @sprintf("%.6f", x)

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
function lineplan_is_feasible(lineplan, network, demand, n_nodes, tp;
                              min_nodes::Int,
                              max_nodes::Int)

    return is_feasible(
        lineplan,
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

# Summary statistics helper
function safe_stats(x::Vector{Float64}; tol::Float64 = 1e-6)
    isempty(x) && return (min = missing, max = missing, avg = missing, median = missing)

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
        o == d && continue

        if isfinite(shortest[o, d])
            total_served_passengers += demand[o, d]
        else
            push!(inf_pairs, (o, d))
            demand[o, d] > 0 && push!(inf_pairs_with_demand, (o, d))
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

# Helper to push proposed iteration history entries
function push_proposed_history!(history,
                                iter,
                                phase,
                                k,
                                move,
                                move_name,
                                method_used,
                                base_lineplan_before,
                                cand_lineplan,
                                base_obj_before,
                                best_obj_before,
                                cand_obj,
                                accepted_flag,
                                improved_flag,
                                current_updated_flag,
                                rejected_flag,
                                used_fallback,
                                detection_time)

    history === nothing && return

    push!(history, (
        iter = iter,
        phase = phase,
        shaking_neighbourhood = k,
        move = deepcopy(move),
        move_type = move_name,
        eval_method_used = method_used,
        current_lineplan_before = deepcopy(base_lineplan_before),
        lineplan = deepcopy(cand_lineplan),
        current_obj_before = base_obj_before,
        local_obj = cand_obj,
        candidate_obj = cand_obj,
        delta_vs_base = cand_obj - base_obj_before,
        delta_vs_best_before = cand_obj - best_obj_before,
        accepted = accepted_flag,
        improved = improved_flag,
        current_updated = current_updated_flag,
        rejected = rejected_flag,
        used_fallback = used_fallback,
        detection_time = detection_time
    ))

    return
end

# First-improvement intensification
function proposed_first_improvement_local_search(start_sol::TransitSolution,
                                                 best_sol::TransitSolution,
                                                 network,
                                                 demand,
                                                 tp,
                                                 n_nodes::Int,
                                                 rng::AbstractRNG;
                                                 eval_mode::Symbol,
                                                 min_nodes::Int = 2,
                                                 max_nodes::Int = typemax(Int),
                                                 max_tries::Int = 10_000,
                                                 max_passes::Int = 100,
                                                 max_failures::Int = 50,
                                                 full_threshold::Float64 = 0.70,
                                                 history = nothing,
                                                 iter_ref::Base.RefValue{Int} = Ref(0),
                                                 time_start::Float64 = time(),
                                                 time_limit_sec::Float64 = Inf)

    current = deepcopy(start_sol)
    best = deepcopy(best_sol)

    accepted = 0
    improved = 0
    rejected = 0
    feasible_candidates = 0
    detection_time_total = 0.0

    passes = 0
    failures = 0

    while passes < max_passes &&
          failures < max_failures &&
          (time() - time_start) < time_limit_sec

        passes += 1
        found_improvement_in_pass = false
        evaluated_candidate_in_pass = false

        # Sequential intensification neighbourhood order
        for fn in PROPOSED_INTENSIFICATION_FNS

            (time() - time_start) < time_limit_sec || break

            cand_lineplan, move, move_name =
                fn(
                    current.lineplan,
                    network,
                    rng;
                    min_nodes = min_nodes,
                    max_nodes = max_nodes
                )

            cand_lineplan === nothing && continue
            changed_lineplan(current.lineplan, cand_lineplan) || continue

            if !candidate_is_feasible(
                    cand_lineplan,
                    network,
                    demand,
                    tp,
                    n_nodes;
                    min_nodes = min_nodes,
                    max_nodes = max_nodes
               )
                continue
            end

            feasible_candidates += 1
            iter_ref[] += 1
            evaluated_candidate_in_pass = true

            base_lineplan_before = deepcopy(current.lineplan)
            base_obj_before = current.obj
            best_obj_before = best.obj

            cand, used_fallback, detection_time, method_used =
                build_solution_proposed_with_meta(
                    cand_lineplan,
                    current,
                    move,
                    eval_mode,
                    network,
                    demand,
                    tp;
                    full_threshold = full_threshold
                )

            detection_time_total += detection_time

            accepted_flag = false
            improved_flag = false
            current_updated_flag = false
            rejected_flag = false

            if cand.obj < current.obj
                current = deepcopy(cand)

                accepted += 1
                accepted_flag = true
                current_updated_flag = true

                failures = 0
                found_improvement_in_pass = true

                if cand.obj < best.obj
                    best = deepcopy(cand)
                    improved += 1
                    improved_flag = true
                end

                push_proposed_history!(
                    history,
                    iter_ref[],
                    "Intensification",
                    missing,
                    move,
                    move_name,
                    method_used,
                    base_lineplan_before,
                    cand_lineplan,
                    base_obj_before,
                    best_obj_before,
                    cand.obj,
                    accepted_flag,
                    improved_flag,
                    current_updated_flag,
                    rejected_flag,
                    used_fallback,
                    detection_time
                )

                # First improvement: accept immediately and restart from first neighbourhood
                break

            else
                rejected += 1
                rejected_flag = true
                failures += 1

                push_proposed_history!(
                    history,
                    iter_ref[],
                    "Intensification",
                    missing,
                    move,
                    move_name,
                    method_used,
                    base_lineplan_before,
                    cand_lineplan,
                    base_obj_before,
                    best_obj_before,
                    cand.obj,
                    accepted_flag,
                    improved_flag,
                    current_updated_flag,
                    rejected_flag,
                    used_fallback,
                    detection_time
                )

                # Important:
                # Do not stop the pass. Try the next intensification move type.
                continue
            end
        end

        if found_improvement_in_pass
            continue
        end

        if !evaluated_candidate_in_pass
            failures += 1
        end
    end

    return (
        current = current,
        best = best,
        accepted = accepted,
        improved = improved,
        rejected = rejected,
        feasible_candidates = feasible_candidates,
        detection_time_total = detection_time_total,
        passes = passes
    )
end

# Main proposed basic VNS function
function proposed_framework(initial_lineplan,
                            network,
                            demand,
                            tp;
                            eval_mode::Symbol = :full,
                            rng::AbstractRNG = MersenneTwister(1),
                            min_nodes::Int = 2,
                            max_nodes::Int = size(initial_lineplan, 2),
                            time_limit_sec::Float64 = 120.0,
                            max_tries::Int = 10_000,
                            shake_max_tries::Int = 20,
                            local_max_passes::Int = 100,
                            local_max_failures::Int = 50,
                            full_threshold::Float64 = 0.70,
                            store_history::Bool = false)

    n_nodes = size(network, 1)

    if !lineplan_is_feasible(
            initial_lineplan,
            network,
            demand,
            n_nodes,
            tp;
            min_nodes = min_nodes,
            max_nodes = max_nodes
       )
        error("Initial lineplan is infeasible.")
    end

    current = build_solution_feasible(initial_lineplan, network, demand, tp)
    best = deepcopy(current)

    iter_ref = Ref(0)

    feasible_candidates = 0
    accepted = 0
    improved = 0
    rejected = 0

    detection_time_total = 0.0
    history = store_history ? NamedTuple[] : nothing

    start_time = time()
    K = proposed_shaking_neighbourhood_count()

    while (time() - start_time) < time_limit_sec
        k = 1

        while k <= K && (time() - start_time) < time_limit_sec

            base_lineplan_before = deepcopy(current.lineplan)
            base_obj_before = current.obj
            best_obj_before = best.obj

            shake_lineplan, shake_move, shake_move_name =
                generate_feasible_proposed_shake_neighbor(
                    current,
                    k,
                    network,
                    demand,
                    tp,
                    n_nodes,
                    rng;
                    min_nodes = min_nodes,
                    max_nodes = max_nodes,
                    max_tries = shake_max_tries
                )

            if shake_lineplan === nothing
                k += 1
                continue
            end

            feasible_candidates += 1
            iter_ref[] += 1

            shaken_sol, used_fallback, detection_time, method_used =
                build_solution_proposed_with_meta(
                    shake_lineplan,
                    current,
                    shake_move,
                    eval_mode,
                    network,
                    demand,
                    tp;
                    full_threshold = full_threshold
                )

            detection_time_total += detection_time

            push_proposed_history!(
                history,
                iter_ref[],
                "Shaking",
                k,
                shake_move,
                shake_move_name,
                method_used,
                base_lineplan_before,
                shake_lineplan,
                base_obj_before,
                best_obj_before,
                shaken_sol.obj,
                false,
                false,
                false,
                false,
                used_fallback,
                detection_time
            )

            local_result =
                proposed_first_improvement_local_search(
                    shaken_sol,
                    best,
                    network,
                    demand,
                    tp,
                    n_nodes,
                    rng;
                    eval_mode = eval_mode,
                    min_nodes = min_nodes,
                    max_nodes = max_nodes,
                    max_tries = max_tries,
                    max_passes = local_max_passes,
                    max_failures = local_max_failures,
                    full_threshold = full_threshold,
                    history = history,
                    iter_ref = iter_ref,
                    time_start = start_time,
                    time_limit_sec = time_limit_sec
                )

            local_sol = local_result.current
            local_best = local_result.best

            feasible_candidates += local_result.feasible_candidates
            accepted += local_result.accepted
            improved += local_result.improved
            rejected += local_result.rejected
            detection_time_total += local_result.detection_time_total

            if local_sol.obj < current.obj
                current = deepcopy(local_sol)

                if current.obj < best.obj
                    best = deepcopy(current)
                elseif local_best.obj < best.obj
                    best = deepcopy(local_best)
                end

                k = 1
            else
                k += 1
            end
        end
    end

    if current.obj < best.obj
        best = deepcopy(current)
    end

    runtime = time() - start_time

    return (
        eval_mode = eval_mode,
        runtime = runtime,
        iterations = iter_ref[],
        feasible_candidates = feasible_candidates,
        best_objective = best.obj,
        final_objective = current.obj,
        best_lineplan = deepcopy(best.lineplan),
        final_lineplan = deepcopy(current.lineplan),
        accepted = accepted,
        improved = improved,
        rejected = rejected,
        avg_time_per_iteration = runtime / max(iter_ref[], 1),
        iterations_per_second = iter_ref[] / max(runtime, 1e-12),
        detection_time_total = detection_time_total,
        avg_detection_time = detection_time_total / max(feasible_candidates, 1),
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

    return (
        tested = tested,
        mismatch_count = mismatch_count,
        mismatch_rate_pct = tested == 0 ? 0.0 : clean_zero(100 * mismatch_count / tested; tol = tol)
    )
end

# Iteration accuracy log
function build_accuracy_log_from_history(history,
                                         network,
                                         demand,
                                         tp;
                                         tol::Float64 = 1e-6,
                                         detection_tol::Float64 = 1e-9,
                                         max_rows::Union{Int,Nothing} = nothing)

    total_demand = sum(demand)

    rows = DataFrame(
        Iteration = Union{Missing,Int}[],
        Phase = String[],
        ShakingNeighbourhood = Any[],
        MoveType = String[],
        EvalMethodUsed = String[],

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
        CurrentUpdated = Union{Missing,Bool}[],
        Rejected = Union{Missing,Bool}[],
        UsedFallback = Union{Missing,Bool}[],
        DetectionTimeSec = String[],

        AnyInfOD = Union{Missing,Bool}[],
        AnyInfODWithDemand = Union{Missing,Bool}[],
        NInfODWithDemand = Union{Missing,Int}[],
        InfODWithDemandPairs = String[],
        TotalServedPassengers = String[]
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
        local_avt = isfinite(local_obj) ? clean_zero(local_obj / total_demand; tol = tol) : Inf
        full_avt = isfinite(full_obj_clean) ? clean_zero(full_obj_clean / total_demand; tol = tol) : Inf

        quality_loss = clean_zero(local_obj - full_obj_clean; tol = tol)
        quality_loss_pct =
            abs(full_obj_clean) <= tol ? 0.0 :
            clean_zero(100 * abs(quality_loss) / abs(full_obj_clean); tol = tol)

        push!(quality_loss_vals, quality_loss)
        push!(quality_loss_pct_vals, quality_loss_pct)
        push!(detection_time_vals, clean_zero(entry.detection_time; tol = tol))

        diag = solution_diagnostics(entry.lineplan, network, demand, tp)

        push!(rows, (
            entry.iter,
            entry.phase,
            entry.shaking_neighbourhood,
            entry.move_type,
            String(entry.eval_method_used),

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
            entry.current_updated,
            entry.rejected,
            entry.used_fallback,
            fmt6(clean_zero(entry.detection_time; tol = tol)),

            diag.any_inf_od,
            diag.any_inf_od_with_demand,
            diag.n_inf_od_with_demand,
            diag.inf_od_with_demand_pairs,
            fmt2(clean_zero(diag.total_served_passengers; tol = tol))
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
        missing,
        "Average",
        "",

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
        missing,
        missing,
        missing,
        fmt6_or_missing(avg_detection_time),

        missing,
        missing,
        missing,
        "",
        ""
    ))

    return rows
end

# Run comparison
function run_proposed_framework_comparison(initial_lineplan,
                                           network,
                                           demand,
                                           tp;
                                           eval_modes::Vector{Symbol},
                                           seed::Int = 1,
                                           n_runs::Int = 10,
                                           min_nodes::Int = 2,
                                           max_nodes::Int = size(initial_lineplan, 2),
                                           time_limit_sec::Float64 = 120.0,
                                           max_tries::Int = 10_000,
                                           shake_max_tries::Int = 20,
                                           local_max_passes::Int = 100,
                                           local_max_failures::Int = 50,
                                           full_threshold::Float64 = 0.70,
                                           tol::Float64 = 1e-6,
                                           build_iteration_logs::Bool = false,
                                           detection_tol::Float64 = 1e-9,
                                           max_log_rows::Int = 10000)

    n_nodes = size(network, 1)
    total_demand = sum(demand)

    if !lineplan_is_feasible(
            initial_lineplan,
            network,
            demand,
            n_nodes,
            tp;
            min_nodes = min_nodes,
            max_nodes = max_nodes
       )
        error("Initial lineplan passed to run_proposed_framework_comparison is infeasible.")
    end

    initial_sol = build_solution_feasible(initial_lineplan, network, demand, tp)
    initial_diag = solution_diagnostics(initial_lineplan, network, demand, tp)

    vns_rows = DataFrame(
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

        FeasibleCandidates_Min = String[],
        FeasibleCandidates_Max = String[],
        FeasibleCandidates_Avg = String[],
        FeasibleCandidates_Median = String[],

        AcceptedMoves_Min = String[],
        AcceptedMoves_Max = String[],
        AcceptedMoves_Avg = String[],
        AcceptedMoves_Median = String[],

        ImprovingMoves_Min = String[],
        ImprovingMoves_Max = String[],
        ImprovingMoves_Avg = String[],
        ImprovingMoves_Median = String[],

        RejectedMoves_Min = String[],
        RejectedMoves_Max = String[],
        RejectedMoves_Avg = String[],
        RejectedMoves_Median = String[]
    )

    vns_run_rows = DataFrame(
        Method = String[],
        Run = Int[],

        RuntimeSec = String[],
        Iterations = Int[],
        FeasibleCandidates = Int[],

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
        ImprovingMoves = Int[],
        RejectedMoves = Int[]
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
    println("PROPOSED BASIC VNS FRAMEWORK COMPARISON")
    @printf("Initial TTT                  = %.6f\n", initial_sol.obj)
    @printf("Initial AVT                  = %.6f\n", initial_sol.obj / total_demand)
    @printf("Initial total served demand  = %.2f\n", initial_diag.total_served_passengers)
    println("Initial AnyInfOD             = ", initial_diag.any_inf_od)
    println("Initial AnyInfODWithDemand   = ", initial_diag.any_inf_od_with_demand)
    println("Initial NInfODWithDemand     = ", initial_diag.n_inf_od_with_demand)
    @printf("Time limit per run           = %.2f s\n", time_limit_sec)
    println("Number of reported runs      = ", n_runs)
    println("Methods                      = ", eval_modes)

    for mode in eval_modes
        println("--------------------------------------------------")
        println("Running proposed Basic VNS for method: ", mode)

        println("   Warm-up run (not reported)")

        proposed_framework(
            initial_lineplan,
            network,
            demand,
            tp;
            eval_mode = mode,
            rng = MersenneTwister(seed),
            min_nodes = min_nodes,
            max_nodes = max_nodes,
            time_limit_sec = time_limit_sec,
            max_tries = max_tries,
            shake_max_tries = shake_max_tries,
            local_max_passes = local_max_passes,
            local_max_failures = local_max_failures,
            full_threshold = full_threshold,
            store_history = false
        )

        runtimes = Float64[]
        iterations_list = Float64[]
        feasible_candidates_list = Float64[]

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
        rejected_list = Float64[]

        best_run_obj = Inf
        best_run_lineplan = deepcopy(initial_lineplan)

        total_tested = 0
        total_mismatch = 0

        for run_id in 1:n_runs
            run_seed = seed + run_id

            println("   Run ", run_id, "/", n_runs, " (seed = ", run_seed, ")")

            result = proposed_framework(
                initial_lineplan,
                network,
                demand,
                tp;
                eval_mode = mode,
                rng = MersenneTwister(run_seed),
                min_nodes = min_nodes,
                max_nodes = max_nodes,
                time_limit_sec = time_limit_sec,
                max_tries = max_tries,
                shake_max_tries = shake_max_tries,
                local_max_passes = local_max_passes,
                local_max_failures = local_max_failures,
                full_threshold = full_threshold,
                store_history = build_iteration_logs
            )

            reported_best_obj = clean_zero(result.best_objective; tol = tol)
            reported_final_obj = clean_zero(result.final_objective; tol = tol)

            best_full_obj, _, _, _, _ =
                evaluate_full_solution(network, result.best_lineplan, demand, tp)

            final_full_obj, _, _, _, _ =
                evaluate_full_solution(network, result.final_lineplan, demand, tp)

            best_full_obj = clean_zero(best_full_obj; tol = tol)
            final_full_obj = clean_zero(final_full_obj; tol = tol)

            best_gap = clean_zero(reported_best_obj - best_full_obj; tol = tol)
            final_gap = clean_zero(reported_final_obj - final_full_obj; tol = tol)

            best_gap_pct =
                abs(best_full_obj) <= tol ? 0.0 :
                clean_zero(100 * abs(best_gap) / abs(best_full_obj); tol = tol)

            final_gap_pct =
                abs(final_full_obj) <= tol ? 0.0 :
                clean_zero(100 * abs(final_gap) / abs(final_full_obj); tol = tol)

            reported_best_avt = clean_zero(reported_best_obj / total_demand; tol = tol)
            full_best_avt = clean_zero(best_full_obj / total_demand; tol = tol)
            reported_final_avt = clean_zero(reported_final_obj / total_demand; tol = tol)
            full_final_avt = clean_zero(final_full_obj / total_demand; tol = tol)

            push!(runtimes, result.runtime)
            push!(iterations_list, Float64(result.iterations))
            push!(feasible_candidates_list, Float64(result.feasible_candidates))

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

            push!(accepted_list, Float64(result.accepted))
            push!(improved_list, Float64(result.improved))
            push!(rejected_list, Float64(result.rejected))

            push!(vns_run_rows, (
                String(mode),
                run_id,

                fmt4(clean_zero(result.runtime; tol = tol)),
                result.iterations,
                result.feasible_candidates,

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

                result.accepted,
                result.improved,
                result.rejected
            ))

            if best_full_obj < best_run_obj
                best_run_obj = best_full_obj
                best_run_lineplan = deepcopy(result.best_lineplan)
            end

            if result.history !== nothing
                println("Checking objective accuracy offline for method: ", mode, " | run: ", run_id)

                acc_result = evaluate_accuracy_from_history(
                    result.history,
                    network,
                    demand,
                    tp;
                    tol = tol
                )

                total_tested += acc_result.tested
                total_mismatch += acc_result.mismatch_count
            end

            if build_iteration_logs && run_id == 1 && result.history !== nothing
                iteration_logs[String(mode)] =
                    build_accuracy_log_from_history(
                        result.history,
                        network,
                        demand,
                        tp;
                        tol = tol,
                        detection_tol = detection_tol,
                        max_rows = max_log_rows
                    )
            end
        end

        best_lineplans[String(mode)] = deepcopy(best_run_lineplan)

        runtime_stats = safe_stats(runtimes; tol = tol)
        iter_stats = safe_stats(iterations_list; tol = tol)
        feasible_candidates_stats = safe_stats(feasible_candidates_list; tol = tol)

        reported_best_ttt_stats = safe_stats(reported_best_ttt_list; tol = tol)
        full_best_ttt_stats = safe_stats(full_best_ttt_list; tol = tol)
        best_gap_stats = safe_stats(best_gap_list; tol = tol)
        best_gap_pct_stats = safe_stats(best_gap_pct_list; tol = tol)

        reported_final_ttt_stats = safe_stats(reported_final_ttt_list; tol = tol)
        full_final_ttt_stats = safe_stats(full_final_ttt_list; tol = tol)
        final_gap_stats = safe_stats(final_gap_list; tol = tol)
        final_gap_pct_stats = safe_stats(final_gap_pct_list; tol = tol)

        reported_best_avt_stats = safe_stats(reported_best_avt_list; tol = tol)
        full_best_avt_stats = safe_stats(full_best_avt_list; tol = tol)
        reported_final_avt_stats = safe_stats(reported_final_avt_list; tol = tol)
        full_final_avt_stats = safe_stats(full_final_avt_list; tol = tol)

        accepted_stats = safe_stats(accepted_list; tol = tol)
        improved_stats = safe_stats(improved_list; tol = tol)
        rejected_stats = safe_stats(rejected_list; tol = tol)

        push!(vns_rows, (
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

            fmt2_or_missing(feasible_candidates_stats.min),
            fmt2_or_missing(feasible_candidates_stats.max),
            fmt2_or_missing(feasible_candidates_stats.avg),
            fmt2_or_missing(feasible_candidates_stats.median),

            fmt2_or_missing(accepted_stats.min),
            fmt2_or_missing(accepted_stats.max),
            fmt2_or_missing(accepted_stats.avg),
            fmt2_or_missing(accepted_stats.median),

            fmt2_or_missing(improved_stats.min),
            fmt2_or_missing(improved_stats.max),
            fmt2_or_missing(improved_stats.avg),
            fmt2_or_missing(improved_stats.median),

            fmt2_or_missing(rejected_stats.min),
            fmt2_or_missing(rejected_stats.max),
            fmt2_or_missing(rejected_stats.avg),
            fmt2_or_missing(rejected_stats.median)
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
        initial_total_served_passengers =
            fmt2(clean_zero(initial_diag.total_served_passengers; tol = tol)),
        initialization_time = nothing,
        initial_lineplan = deepcopy(initial_lineplan),
        best_lineplans = best_lineplans,
        vns_results = vns_rows,
        vns_run_results = vns_run_rows,
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
            candidate =
                (length(base_name) > trim_len ? base_name[1:trim_len] : base_name) * suffix

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
                "Initialization time (sec)"
            ],
            Value = [
                results.initial_ttt,
                results.initial_avt,
                string(results.initial_any_inf_od),
                string(results.initial_any_inf_od_with_demand),
                string(results.initial_n_inf_od_with_demand),
                results.initial_total_served_passengers,
                results.initialization_time
            ]
        )

        sheet_summary = XLSX.addsheet!(xf, make_unique_sheet_name("Summary", used_sheet_names))
        XLSX.writetable!(sheet_summary, Tables.columntable(summary_df); write_columnnames = true)

        sheet_results = XLSX.addsheet!(xf, make_unique_sheet_name("VNS_Results", used_sheet_names))
        XLSX.writetable!(sheet_results, Tables.columntable(results.vns_results); write_columnnames = true)

        sheet_runs = XLSX.addsheet!(xf, make_unique_sheet_name("VNS_Run_Results", used_sheet_names))
        XLSX.writetable!(sheet_runs, Tables.columntable(results.vns_run_results); write_columnnames = true)

        sheet_accuracy = XLSX.addsheet!(xf, make_unique_sheet_name("Accuracy_vs_Full", used_sheet_names))
        XLSX.writetable!(sheet_accuracy, Tables.columntable(results.accuracy_results); write_columnnames = true)

        sheet_init = XLSX.addsheet!(xf, make_unique_sheet_name("Initial_Lineplan", used_sheet_names))
        XLSX.writetable!(
            sheet_init,
            Tables.columntable(DataFrame(results.initial_lineplan, :auto));
            write_columnnames = false
        )

        for (method, lp) in results.best_lineplans
            sheet = XLSX.addsheet!(xf, make_unique_sheet_name("Best_" * method, used_sheet_names))
            XLSX.writetable!(
                sheet,
                Tables.columntable(DataFrame(lp, :auto));
                write_columnnames = false
            )
        end

        for (method, df) in results.iteration_logs
            sheet = XLSX.addsheet!(xf, make_unique_sheet_name(replace(method, ":" => ""), used_sheet_names))
            XLSX.writetable!(sheet, Tables.columntable(df); write_columnnames = true)
        end
    end

    println("==================================================")
    println("Results exported to: ", filename)
end

# Data loading
BASE = normpath(joinpath(@__DIR__, ".."))

network_path = joinpath(BASE, "TestNetwork.txt")
demand_path = joinpath(BASE, "TestDemand.txt")
lineplan_path = joinpath(BASE, "TestLineplan.txt")

network = Float64.(readdlm(network_path))
demand = Float64.(readdlm(demand_path))
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

loaded_feasible = lineplan_is_feasible(
    loaded_lineplan,
    network,
    demand,
    n_nodes,
    tp;
    min_nodes = min_nodes,
    max_nodes = max_nodes
)

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

initialized_feasible = lineplan_is_feasible(
    initial_lineplan,
    network,
    demand,
    n_nodes,
    tp;
    min_nodes = min_nodes,
    max_nodes = max_nodes
)

println("Constructed initial lineplan feasible? ", initialized_feasible)

if !initialized_feasible
    error("Constructed initial lineplan is infeasible.")
end

@printf("Initialization time = %.4f seconds\n", init_runtime)

# Methods to compare
eval_modes = [
    :full,
    :local_dijkstra_pq,
    :local_floyd_subset,
    :threshold_basic,
    :adaptive_local
]

# Run experiment
results = run_proposed_framework_comparison(
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
    max_tries = 10_000,
    shake_max_tries = 20,
    local_max_passes = 100,
    local_max_failures = 50,
    full_threshold = 0.70,
    tol = 1e-6,
    build_iteration_logs = false,
    detection_tol = 1e-9,
    max_log_rows = 1000
)

results = merge(results, (initialization_time = fmt4(init_runtime),))

export_comparison_to_excel(
    joinpath(BASE, "proposed_basic_vns_framework_comparison.xlsx"),
    results
)

nothing