using DelimitedFiles, XLSX, DataFrames, Tables

include(joinpath(@__DIR__, "feasibility.jl"))
include(joinpath(@__DIR__, "solution.jl"))
include(joinpath(@__DIR__, "evaluation.jl"))

BASE = normpath(joinpath(@__DIR__, ".."))

network_path  = joinpath(BASE, "TestNetwork.txt")
demand_path   = joinpath(BASE, "TestDemand.txt")
lineplan_path = joinpath(BASE, "TestLineplan.txt")
output_path   = joinpath(BASE, "lineplan_check_results.xlsx")

network  = Float64.(readdlm(network_path))
demand   = Float64.(readdlm(demand_path))
lineplan = readdlm(lineplan_path, Int)

tp = 5.0
min_nodes = 2
max_nodes = size(lineplan, 2)
n_nodes = size(network, 1)

# Helper to convert sets of OD pairs to sorted string representation for diagnostics
function pairset_to_string(pairs)
    isempty(pairs) && return ""
    vecpairs = sort!(collect(pairs), by = x -> (x[1], x[2]))
    return join(["($(o),$(d))" for (o, d) in vecpairs], ", ")
end

# Compute present and absent stops for diagnostics
function compute_present_absent_stops(lineplan, n_nodes)
    present = nodes_present(lineplan, n_nodes)
    absent = .!present
    return present, absent
end

# Check structural feasibility of a lineplan and report the first detected violation
function explain_structure_only(lineplan, network;
                                min_nodes::Int = 2,
                                max_nodes::Int = typemax(Int))

    for L in 1:size(lineplan, 1)
        route = route_from_lineplan(lineplan, L)

        # allow all-zero rows
        if isempty(route)
            continue
        end

        if length(route) < min_nodes || length(route) > max_nodes
            return false, "Line $L has length $(length(route)), which is outside [$min_nodes, $max_nodes]"
        end

        if !is_valid_route(route, network;
                           min_nodes = min_nodes,
                           max_nodes = max_nodes,
                           allow_loop_lines = true)

            if any(!has_edge(network, route[i], route[i + 1]) for i in 1:(length(route) - 1))
                for i in 1:(length(route) - 1)
                    u = route[i]
                    v = route[i + 1]
                    if !has_edge(network, u, v)
                        return false, "Line $L has no network edge between $u and $v"
                    end
                end
            end

            if length(unique(route)) != length(route)
                is_closed_loop = is_proper_loop(route)
                if !is_closed_loop
                    return false, "Line $L contains repeated nodes that do not form a simple loop"
                end
            end

            return false, "Line $L is structurally infeasible"
        end
    end

    return true, "Lineplan is structurally feasible"
end

# Compute objective values and diagnostics for a lineplan, assuming it is structurally feasible
function compute_objective_structure_only(lineplan,
                                          network,
                                          demand,
                                          tp;
                                          min_nodes::Int = 2,
                                          max_nodes::Int = size(lineplan, 2))

    structurally_feasible, _ = explain_structure_only(
        lineplan,
        network;
        min_nodes = min_nodes,
        max_nodes = max_nodes
    )

    shortest, lineinfo, directlink, path = shortest_paths(network, lineplan, tp)
    served_stops = get_served_stops(shortest)

    if !structurally_feasible
        TTT = Inf
        AVT = Inf
    else
        TTT = sum(shortest[served_stops, served_stops] .* demand[served_stops, served_stops])
        AVT = TTT / sum(demand)
    end

    return TTT, AVT, served_stops, shortest, lineinfo, directlink, path
end

# Compute solution diagnostics for OD connectivity, served demand, and stop coverage
function solution_diagnostics(lineplan, network, demand, tp)
    shortest, _, _, _ = shortest_paths(network, lineplan, tp)
    n_nodes = size(network, 1)

    inf_pairs = Set{Tuple{Int,Int}}()
    inf_pairs_with_demand = Set{Tuple{Int,Int}}()

    total_served_passengers = 0.0
    total_od_pairs_with_demand = 0

    for o in 1:n_nodes, d in 1:n_nodes
        if o == d
            continue
        end

        if demand[o, d] > 0
            total_od_pairs_with_demand += 1
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

    present_stops, absent_stops = compute_present_absent_stops(lineplan, n_nodes)

    return (
        any_inf_od = !isempty(inf_pairs),
        any_inf_od_with_demand = !isempty(inf_pairs_with_demand),
        total_od_pairs_with_demand = total_od_pairs_with_demand,
        n_inf_od_with_demand = length(inf_pairs_with_demand),
        inf_od_with_demand_pairs = pairset_to_string(inf_pairs_with_demand),
        total_served_passengers = total_served_passengers,
        n_present_stops = count(present_stops),
        n_absent_stops = count(absent_stops),
        absent_stop_indices = findall(absent_stops)
    )
end

# Run
struct_ok, reason = explain_structure_only(
    lineplan,
    network;
    min_nodes = min_nodes,
    max_nodes = max_nodes
)

println("==========================================")
println("LINEPLAN FEASIBILITY CHECK")
println("Structurally feasible? ", struct_ok)
println("Reason: ", reason)

TTT = missing
AVT = missing
diag = nothing

if struct_ok
    TTT, AVT, served_stops, shortest, lineinfo, directlink, path =
        compute_objective_structure_only(
            lineplan,
            network,
            demand,
            tp;
            min_nodes = min_nodes,
            max_nodes = max_nodes
        )

    diag = solution_diagnostics(lineplan, network, demand, tp)

    println()
    println("OBJECTIVE VALUES")
    println("TTT: ", TTT)
    println("AVT: ", AVT)

    println()
    println("DIAGNOSTICS")
    println("AnyInfOD: ", diag.any_inf_od)
    println("AnyInfODWithDemand: ", diag.any_inf_od_with_demand)
    println("TotalODPairsWithDemand: ", diag.total_od_pairs_with_demand)
    println("NInfODWithDemand: ", diag.n_inf_od_with_demand)
    println("InfODWithDemandPairs: ", diag.inf_od_with_demand_pairs)
    println("TotalServedPassengers: ", diag.total_served_passengers)
    println("NPresentStops: ", diag.n_present_stops)
    println("NAbsentStops: ", diag.n_absent_stops)
    println("AbsentStops: ", diag.absent_stop_indices)
end

# Export to Excel
results_df = DataFrame(
    Metric = String[],
    Value = Any[]
)

push!(results_df, ("Structurally feasible", struct_ok))
push!(results_df, ("Reason", reason))

if struct_ok
    push!(results_df, ("TTT", TTT))
    push!(results_df, ("AVT", AVT))
    push!(results_df, ("AnyInfOD", diag.any_inf_od))
    push!(results_df, ("AnyInfODWithDemand", diag.any_inf_od_with_demand))
    push!(results_df, ("TotalODPairsWithDemand", diag.total_od_pairs_with_demand))
    push!(results_df, ("NInfODWithDemand", diag.n_inf_od_with_demand))
    push!(results_df, ("InfODWithDemandPairs", diag.inf_od_with_demand_pairs))
    push!(results_df, ("TotalServedPassengers", diag.total_served_passengers))
    push!(results_df, ("NPresentStops", diag.n_present_stops))
    push!(results_df, ("NAbsentStops", diag.n_absent_stops))
    push!(results_df, ("AbsentStops", string(diag.absent_stop_indices)))
end

XLSX.openxlsx(output_path, mode = "w") do xf
    sheet = XLSX.addsheet!(xf, "Results")
    XLSX.writetable!(sheet, Tables.columntable(results_df); write_columnnames = true)
end

println()
println("Results exported to: ", output_path)