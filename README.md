# Local Evaluations of Passenger Assignment to Speed up Line Planning Algorithms

In this folder, 

The network instances used in this study are provided in the `NetworkInstances` folder. For each instance, the demand matrix, line plan, and infrastructure network are given in separate `.txt` files. Coordinates for Mandl's Swiss network and the Mumford instances are also included. These instance files were obtained from: https://www.mech.kuleuven.be/en/mim/lp#section-28. 

The implementation files used for local move definitions, impact analysis, partial PAP evaluation, and metaheuristic frameworks are included in the `Codes` folder. 

File organization is as follows:

- Core implementation files:  
  `BasicElementaryMoves.jl`, `LocalMoves_BEM.jl`, `evaluator.jl`

- Impact-analysis experiment files:  
  `RandomExperiments.jl`, `main_LocalMoves.jl`

- Partial PAP evaluation files:  
  `Rules.jl`, `Rules_StraightforwardApproach.jl`, `LocalEvaluation.jl`, sweep files

- Metaheuristic utility files: 
  `common_lineplan_utils.jl`, `initialization.jl`, `feasibility.jl`, `evaluation.jl`, `solution.jl`

- Simulated Annealing (SA) framework:  
  `sa_evaluation.jl`, `sa_moves.jl`, `sa_neighborhood.jl`, `SA_main.jl`

- Variable Neighborhood Search (VNS) framework:  
  `vns_evaluation.jl`, `vns_moves.jl`, `vns_neighborhood.jl`, `VNS_main.jl`

- Proposed basic VNS framework:
  `basic_vns_evaluation.jl`, `basic_vns_moves.jl`, `basic_vns_neighborhood.jl`, `basic_vns_main.jl`

- Analysis tool:  
  `lineplan_check.jl`


## Core implementation files

### `BasicElementaryMoves.jl`:
Contains the functions for the basic elementary moves. 

Included basic elementary moves: 
- `Insert(line, node, position)`: Inserts a single node into a line before a given position,
- `Remove(line, node)`: Removes a single node from a line, 
- `InsertSegment(line, segment, after node)`: Inserts an ordered subsequence of nodes into a line after a given node, 
- `RemoveSegment(line, start node, end node)`: Removes an ordered subsequence of nodes from a line between two given nodes, while retaining the boundary nodes, 
- `Reverse(segment)`: Reverses the order of nodes within a selected segment, 
- `CreateSegment([node,...,node])`: Creates a new segment from an ordered subsequence of nodes.

### `LocalMoves_BEM.jl`:
Contains the functions for the 14 local moves used in the thesis, based on the basic elementary moves defined in `BasicElementaryMoves.jl`. In addition, move-specific feasibility and continuity checks are included.

Included local moves which return the new lineplan obtained after the application of the local move:
- Insertion at Terminal
- Insertion
- Removal at Terminal
- Removal
- Swap within Line
- Swap between Lines
- Position change within line (Transfer within line)
- Position change between lines (Transfer between lines)
- Substitution
- Partial Segment Swap
- Reversal
- Segment Reversal
- Line Addition
- Line Removal

Included feasibility helper functions:
- `is_feasible_lineplan(lineplan, infnetwork)`: checks infrastructure connectivity.
- `serves_all_demand(lineplan, network, tp, demand)`: checks the serviceability of OD pairs with nonzero demand.
- `is_feasible_candidate(lineplan, network, demand, tp)`: performs overall feasibility checking.

### `evaluator.jl`:
Contains the core evaluation functions. This file is based on the IP1 Assignment 1 code used in the academic year 2025–2026. 

The `directlinknetworkinit` function was modified to compute backward-direction ride distances for non-loop lines. Instead of assigning the forward ride distance to both directions, the reverse-direction ride distance is calculated from the opposite ordering of nodes along the line.

Included functions:
- `directlinknetworkinit(network, lineplan)`: given the infrastructure network and the lineplan, constructs the DLN and associated line information.
- `floyd_warshall(Dist)`: computes all-pairs shortest paths using the Floyd–Warshall algorithm given a distance matrix.
- `shortest_paths(network, lineplan, tp)`: evaluates a line plan by building the DLN, incorporating transfer penalties, and computing shortest travel times.


## Impact-analysis experiment files

### `RandomExperiments.jl`:
Contains the experimental framework used for the *Impact Analysis of Local Moves*. Generates random moves under feasibility constraints, evaluates OD impacts, and performs statistical result analysis.

Main components:
1. OD-Change Assessment
    - `assess_od_changes_updated`: compares the old and new shortest-path results for all OD pairs and quantifies changes.
    - `evaluate_move`: evaluates a single local move by comparing the original and modified line plans.

2. Feasibility
    - `did_change`: checks whether a candidate line plan differs from the original line plan.
    - `enforce_feasible_candidate`: accepts a candidate only if it changes the solution and satisfies feasibility conditions. 
    - `enforce_partial_segment_swap_feasible`: Specialized feasibility checker for the partial segment swap move: no repeated nodes within a line.

3. Random Move Generation and Experimentation for Impact Analysis

    For each of the 14 local moves, the following are included:
    - `all_valid_<move>`: enumerates all valid move instances under the current line plan.
    - `random_<move>`: randomly selects one valid move instance.
    - `run_random_experiments_<move>`: performs repeated randomized experiments, evaluates feasible candidates, and records move impacts.

### `main_LocalMoves.jl`:
Calls the functions from `RandomExperiments.jl` to perform randomized experiments for all local moves, evaluate their impact, and export the results to the Excel file `LocalMoves_Results.xlsx`.


## Partial PAP evaluation files

### `Rules.jl`:
Contains move-specific detection rules for identifying the potentially affected direct links and the functions implemented to derive the potentially affected OD pairs based on the set of detected direct links.

Main helper functions:
- `dynamic_improvement_pairs_tp`: detects OD pairs affected by improved or newly created direct links.
- `dynamic_witness_pairs_tp`: detects OD pairs affected by worsened or removed direct links, or by line-assignment changes.

Included move-specific detection rules:
- `detect_insertion_terminal`
- `detect_insertion_internal`
- `detect_removal_terminal`
- `detect_removal_internal`
- `detect_swap_within_line`
- `detect_swap_between_lines`
- `detect_transfer_within_line`
- `detect_transfer_between_lines`
- `detect_substitution`
- `detect_partial_segment_swap`
- `detect_reversal`
- `detect_segment_reversal`
- `detect_line_addition`
- `detect_line_removal`

Each detection rule compares the old and new direct-link networks after a local move, identifies move-specific direct-link changes based on the associated detection rules, and returns the set of potentially affected OD pairs using the OD-pair identification functions.

### `Rules_StraightforwardApproach.jl`:
Contains move-specific detection rules for identifying affected OD pairs based on direct-link changes for the *Straightforward Approach*. 

This file is mainly an adaptation of `Rules.jl`. The main difference lies in the detection procedure used to identify potentially affected OD pairs from the set of detected direct links.

Main helper function:
- `direct_link_endpoint_pairs_symmetric(edges)`: converts detected direct-link changes into endpoint-based affected OD pairs.

The same move-specific detection rules are included. The only difference is that the detected set of direct links is forwarded to the endpoint-based derivation procedure.

Each detection rule compares the old and new direct-link networks after a local move, identifies move-specific direct-link changes based on the associated detection rules, and returns the set of potentially affected OD pairs using the OD-pair identification functions.

### `LocalEvaluation.jl`:
Contains the PAP evaluation methods used to assess local evaluation approaches and compare them against full evaluation for the *Performance of Partial PAP*.

Main ground-truth reconstruction function:
- `affected_od_groundtruth`: computes the ground-truth affected OD pairs based on changes in travel cost, stop sequence, or line sequence.

Included helper functions:
- `add_transfer_penalty(directlink, tp)`: adds transfer penalties to finite direct-link costs.
- `subtract_first_boarding!(shortest, tp)`: removes the first-boarding transfer penalty from shortest-path matrices.
- `group_by_origin(pairs)`: groups OD pairs by origin node.
- `costmatrix_to_edge_list(cost)`: converts a cost matrix into an edge-list representation.
- `nodes_from_od_pairs(A_rule)`: extracts nodes appearing in detected OD pairs.
- `subset_nodes_from_affected_local`: builds the node subset used for subset-based Floyd–Warshall evaluation.

Included evaluation methods:
- `full_evaluation_split`: performs full evaluation (DLN construction + Floyd–Warshall) with separate timing measurements.
- `full_evaluation`: performs full evaluation and returns total computation time.
- `local_dijkstra_evaluation`: performs local evaluation using basic Dijkstra.
- `local_dijkstra_pq_evaluation`: performs local evaluation using priority-queue Dijkstra (binary-heap).
- `local_fw_subset_evaluation`: performs local evaluation using subset-based Floyd–Warshall.

Included shortest-path algorithms:
- `dijkstra_matrix(cost, source)`: basic Dijkstra implementation.
- `dijkstra_priority_queue(edges_by_origin, source)`: priority-queue Dijkstra (binary-heap) implementation.

Included result analysis functions:
- `detection_quality(A_true, A_rule)`: compares detected affected OD pairs against ground truth.
- `subset_quality`: evaluates subset-based Floyd–Warshall accuracy relative to full evaluation.

### Sweep files:
Contain move-specific sweep experiments used to evaluate local PAP evaluation methods for individual local moves for the *Performance of Partial PAP*.

Each sweep file follows the same general structure:
- resolves the project root and includes the required source files,
- loads the network and demand data, 
- enumerates feasible move candidates for one specific local move,
- computes the full evaluation as ground truth,
- applies the corresponding detection rule,
- compares detected affected OD pairs with the ground truth,
- evaluates local basic Dijkstra, local priority-queue Dijkstra (binary-heap), and local subset-based Floyd–Warshall methods,
- computes objective-value gaps and timing results,
- exports the results to an Excel file.

Corresponding files are:
- `RunInsertion.jl`
- `RunInsertionTerminal.jl`
- `RunLineAddition.jl`
- `RunLineRemoval.jl`
- `RunPartialSegmentSwap.jl`
- `RunRemoval.jl`
- `RunRemovalTerminal.jl`
- `RunReversal.jl`
- `RunSegmentReversal.jl`
- `RunSubstitution.jl`
- `RunSwapBetweenLines.jl`
- `RunSwapWithinLine.jl`
- `RunTransferBetweenLines.jl`
- `RunTransferWithinLine.jl`

The generated Excel output contains detection quality, objective gaps, subset accuracy, and computation-time comparisons for the corresponding local move.


## Metaheuristic utility files

### `common_lineplan_utils.jl`:
Contains helper functions shared across metaheuristic frameworks.

Included helper functions:
- `has_edge(network, u, v)`: checks infrastructure connectivity.
- `route_from_lineplan(lineplan, L)`: extracts the nonzero stop sequence of a line from the line plan.
- `line_stops(lineplan, line_id)`: returns the stop sequence of a line.
- `set_route!(lineplan, line_id, route)`: overwrites a line route in the line plan using a new stop sequence.
- `is_proper_loop(route)`: checks whether a route forms a valid loop line.
- `is_valid_route`: validates a route against node-count, uniqueness/loop, and network-connectivity conditions.

These helper functions are defined conditionally (`if !isdefined(...)`) to avoid redefinition conflicts when shared across multiple files.

### `initialization.jl`:
Contains the initialization framework used to generate feasible initial line plans for the metaheuristic algorithms.

Included helper functions:
- `neighbors_in_network(network, u)`: returns all neighboring nodes connected to node `u`.
- `problematic_od_pairs(lineplan, network, demand, tp)`: identifies OD pairs with positive demand but no feasible path.

Included feasibility function:
- `is_feasible_lineplan`: checks whether a lineplan satisfies feasibility conditions.

Included route construction functions:
- `random_connected_route`: generates a random connected route satisfying route constraints.
- `construct_initial_lineplan`: constructs an initial lineplan using randomized route generation.

Included OD coverage repair functions:
- `try_repair_problematic_od`: attempts to repair infeasible OD coverage by replacing a line route.
- `repair_od_coverage`: iteratively repairs OD coverage violations in a lineplan.

Included initialization functions:
- `initialize_feasible_lineplan`: generates a feasible initial lineplan using construction and repair procedures.
- `initialize_feasible_solution`: generates an initial feasible `TransitSolution` object.

### `feasibility.jl`:
Contains the feasibility functions used to validate line plans in the metaheuristic frameworks.

Included helper functions:
- `covered_nodes(lineplan)`: returns the set of nodes covered by the lineplan.
- `nodes_present(lineplan, n_nodes)`: returns a Boolean vector indicating which nodes appear in the lineplan.
- `demand_nodes(demand, n_nodes)`: returns the set of nodes involved in nonzero demand OD pairs.
- `demand_nodes_present(demand, n_nodes)`: returns a Boolean vector indicating nodes involved in positive-demand OD pairs.
- `build_route_set_graph(lineplan, n_nodes)`: constructs the route-set graph induced by the lineplan.

Included feasibility functions:
- `is_connected_route_set(lineplan, demand, n_nodes)`: checks whether the route-set graph is connected over demand-relevant nodes.
- `serves_all_od_pairs(lineplan, network, demand, n_nodes, tp)`: checks whether all nonzero demand OD pairs are served by the lineplan.
- `is_feasible`: validates a lineplan against route validity, OD coverage, and connectivity requirements.

### `evaluation.jl`:
Contains a unified evaluation framework used for both SA and VNS metaheuristics. 

Included objective and solution helpers:
- `get_served_stops(shortest)`: identifies served stops from a shortest-path matrix.
- `make_solution`: creates a `TransitSolution` object.
- `compute_objective_from_shortest`: computes objective values and path information from shortest-path results.
- `evaluate_full_solution`: performs full solution evaluation.
- `build_solution`: constructs a `TransitSolution` from a lineplan evaluation.

Included detection and preprocessing helpers:
- `detect_affected_ods`: applies the appropriate move-specific detection rule.
- `affected_pairs_by_origin`: groups affected OD pairs by origin.
- `removed_terminal_node`: extracts the removed terminal node for terminal-removal moves.
- `build_directlink_tp`: constructs transfer-penalized direct-link matrices.

Included path reconstruction and shortest-path helpers:
- Floyd–Warshall path encoding and reconstruction utilities (`encode_path_as_fw_split!`, `update_fw_path_for_od!`, `reconstruct_fw_path_nodes`).
- Dijkstra implementations (`dijkstra_basic_tp`, `dijkstra_pq_tp`).
- Subset-based Floyd–Warshall utilities (`floyd_warshall_subset_tp`, `subset_nodes_from_affected`, `induced_submatrix`).

Included local evaluation methods:
- `evaluate_local_dijkstra_*`: local Dijkstra-based evaluation.
- `evaluate_local_floyd_subset_*`: subset-based Floyd–Warshall local evaluation.
- `evaluate_local_hybrid_*`: hybrid local evaluation combining multiple methods.
- `evaluate_local_threshold_adaptive_*`: threshold-adaptive evaluation selecting among local and full evaluation strategies.

Included solution builders:
- `build_solution_local_*`: wrappers constructing evaluated `TransitSolution` objects for individual evaluation modes.
- `build_solution`: unified solution builder for all evaluation modes.
- `build_solution_with_meta`: unified solution builder for all evaluation modes returning additional evaluation metadata.

Evaluation modes include:
- `:full`,
- `:local_dijkstra_basic`,
- `:local_dijkstra_pq`,
- `:local_floyd_subset`,
- `:local_floyd_subset_path`,
- `:local_fw_subset_path_basic`,
- `:local_fw_subset_path_pq`,
- `:local_hybrid_basic`,
- `:local_hybrid_pq`,
- `:local_threshold_basic`,
- `:local_threshold_pq`,
- `:local_floyd_subset_periodic`,
- `:local_hybrid_basic_periodic`,
- `:local_hybrid_pq_periodic`

Additionally, `:adaptive_local` is also included

### `solution.jl`:
Defines the `TransitSolution` structure used to store evaluated candidate solutions across the metaheuristic frameworks.

Stored solution components:
- `lineplan`: transit line plan.
- `obj`: objective value.
- `shortest`: shortest-path travel-time matrix.
- `directlink`: direct-link travel-time matrix.
- `path_old`: stored shortest-path representation used by local evaluation methods.
- `lineinfo`: direct-link line-assignment information.


## Simulated Annealing (SA) framework

### `sa_evaluation.jl`:
Contains the SA evaluation module. This file serves as a wrapper around the unified evaluation framework by importing the shared evaluation functions from `evaluation.jl`. The module enables SA algorithm to use the common evaluation framework.

### `sa_moves.jl`:
Contains the moves definitions and move-application functions used in the SA framework.

Included move types:
- Insertion at terminal
- Removal at terminla
- Reversal

### `sa_neighborhood.jl`:
Contains the neighboring-solution generation procedures used in the SA framework.

Main helper functions:
- `candidate_terminal_insertions`: identifies feasible terminal insertion candidates connected to the last node of a route.

Included neighbor-generation functions:
- `make_small_change`: generates a neighboring line plan by applying a small randomized modification.
- `generate_feasible_neighbor`: repeatedly generates candidate neighbors until a feasible line plan is obtained or the maximum number of attempts is reached.

The file uses helper utilities from `common_lineplan_utils.jl` and applies move operators from `LocalMoves_BEM.jl`.

Generated neighbors are returned for subsequent evaluation and affected OD-pair detection.

### `SA_main.jl`:
Contains the main SA framework. Integrates initialization, neighborhood generation, evaluation, diagnostics, offline accuracy analysis, multi-method comparison, and result export procedures.

Main helper functions:
- `clean_zero`, `fmt*`: numerical cleanup and formatting utilities for reporting and Excel export.
- `lineplan_is_feasible`: wrapper around feasibility checks for SA framework.
- `solution_diagnostics`: analyzes infeasible travel times and served demand.
- `accept_move`: implements the SA acceptance criterion.
- `estimate_sa_schedule_paper`: estimates SA temperature schedule parameters.
- `move_type_name`: converts move objects to readable names.

Included optimization functions:
- `simulated_annealing`: executes the main SA algorithm using the selected evaluation mode and neighborhood generator.
- `run_multi_method_comparison`: runs repeated experiments across multiple evaluation methods and collects performance statistics.

Included accuracy-analysis functions (if set as `true`):
- `evaluate_accuracy_from_history`: performs offline comparison between local and full evaluation objectives.
- `build_accuracy_log_from_history`: constructs detailed iteration-level accuracy logs.

Included export functions:
- `export_comparison_to_excel`: exports summary tables, run statistics, accuracy results, lineplans, and iteration logs to Excel.


## Variable Neighborhood Search (VNS) framework

### `vns_evaluation.jl`:
Contains the VNS evaluation module. This file serves as a wrapper around the unified evaluation framework by importing the shared evaluation functions from `evaluation.jl`. The module enables VNS algorithm to use the common evaluation framework.

### `vns_moves.jl`:
Contains the moves definitions and move-application functions used in the VNS framework.

Included move types:
- Insertion
- Removal 
- Swap within line 
- Substitution 
- Partial line swap 
- Reversal 
- Insertion at terminal 
- Removal at terminal 

### `vns_neighborhood.jl`:
Contains the neighboring-solution generation procedures used in the VNS framework.

Included helper functions:
- `vns_line_stops`: extracts the stop sequence of a line.
- `candidate_nodes_not_in_route`, `unordered_node_pairs`, `consecutive_pairs`: generate candidate nodes and node pairs used for neighborhood construction.
- `candidate_is_feasible`: verifies feasibility of a candidate line plan.
- `safe_insertion_at_terminal`, `safe_removal_at_terminal`: safe wrappers for terminal move operations.

Included neighborhood generation functions:
- `random_vns_neighbor`: randomly selects a neighborhood type, generates a candidate line plan, and validates feasibility.
- `generate_feasible_vns_neighbor`: repeatedly generates candidate neighbors until a feasible line plan is obtained or the maximum number of attempts is reached.

The file uses helper utilities from `common_lineplan_utils.jl` and applies move operators from `LocalMoves_BEM.jl`.

Generated neighbors are returned for subsequent evaluation and affected OD-pair detection.

### `VNS_main.jl`:
Contains the main VNS framework. Integrates initialization, neighborhood generation, evaluation, diagnostics, offline accuracy analysis, multi-method comparison, and result export procedures.

Main helper functions:
- `clean_zero`, `fmt*`: numerical cleanup and formatting utilities for reporting and Excel export.
- `lineplan_is_feasible`: wrapper around feasibility checks for VNS framework.
- `solution_diagnostics`: analyzes infeasible travel times and served demand.
- `use_best_block`: alternates between best and current solutions during neighborhood generation.

Included optimization functions:
- `variable_neighborhood_search`: executes the main VNS algorithm using the selected evaluation mode and neighborhood generator.
- `run_multi_method_comparison`: runs repeated experiments across multiple evaluation methods and collects performance statistics.

Included accuracy-analysis functions (if set as `true`):
- `evaluate_accuracy_from_history`: performs offline comparison between local and full evaluation objectives.
- `build_accuracy_log_from_history`: constructs detailed iteration-level accuracy logs.

Included export functions:
- `export_comparison_to_excel`: exports summary tables, run statistics, accuracy results, lineplans, and iteration logs to Excel.


## Proposed basic VNS framework:

### `basic_vns_evaluation.jl`:
Contains the evaluation framework used in the proposed basic VNS metaheuristic. 

Each detection function dispatches to the corresponding analytical detection rule from `Rules.jl` and returns the set of potentially affected OD pairs.

Included helper functions:
- `line_stops_for_subset`: collects affected line stops before and after a move for subset-based local evaluation.
- `is_loop_line`: checks whether a route contains repeated nodes.

Included evaluation-selection functions:
- `build_solution_proposed_adaptive`: adaptively selects the evaluation strategy depending on move type, loop-line structure, and subset characteristics.
- `build_solution_proposed_with_meta`: general wrapper that dispatches to the selected evaluation mode.

Evaluation modes include:
- `:full`
- `:local_dijkstra_pq`
- `:local_floyd_subset`
- `:threshold_basic`
- `:adaptive_local`

### `basic_vns_moves.jl`:
Contains the moves definitions and move-application functions used in the proposed basic VNS framework.

Included move types:
- Insertion
- Removal
- Swap within line
- Substitution
- Reversal
- Insertion at terminal
- Removal at terminal
- Segment reversal
- Swap between lines
- Position change between lines (Transfer between lines)
- Partial line swap

### `basic_vns_neighborhood.jl`:
Contains the neighboring-solution generation procedures used in the proposed basic VNS framework.

Main helper functions:
- `proposed_line_stops`: extracts the stop sequence of a line.
- `candidate_nodes_not_in_route`, `unordered_node_pairs`, `consecutive_pairs`: generate candidate nodes and node pairs used for neighborhood construction.
- `candidate_is_feasible`: verifies feasibility of a candidate line plan.
- `safe_insertion_at_terminal`, `safe_removal_at_terminal`: safe wrappers for terminal move operations.

Included neighbor-generation functions:

- `proposed_random_intensification_neighbor`: generates an intensification neighbor using the proposed neighborhood operators.
- `proposed_shaking_neighbourhood_count`: returns the number of available shaking neighborhoods.
- `proposed_random_shake_neighbor`: generates a shaking neighbor from the selected neighborhood structure.
- `generate_feasible_proposed_intensification_neighbor`: repeatedly generates intensification neighbors until a feasible candidate is obtained or the maximum number of attempts is reached.
- `generate_feasible_proposed_shake_neighbor`: repeatedly generates shaking neighbors until a feasible candidate is obtained or the maximum number of attempts is reached.

The file uses helper utilities from `common_lineplan_utils.jl` and applies move operators from `LocalMoves_BEM.jl`.

Generated neighbors are returned for subsequent evaluation and affected OD-pair detection.

### `basic_vns_main.jl`:
Contains the main proposed basic VNS framework. Integrates initialization, neighborhood generation, evaluation, diagnostics, offline accuracy analysis, multi-method comparison, and result export procedures.

Main helper functions:
- `clean_zero`, `fmt*`: numerical cleanup and formatting utilities for reporting and Excel export.
- `lineplan_is_feasible`: wrapper around feasibility checks for proposed basic VNS framework.
- `solution_diagnostics`: analyzes infeasible travel times and served demand.
- `push_proposed_history!`: stores iteration-level diagnostics and evaluation metadata.

Included optimization functions:
- `proposed_first_improvement_local_search`: performs first-improvement intensification using sequential neighbourhood exploration.
- `proposed_framework`: executes the main proposed basic VNS algorithm using shaking, local search, and the selected evaluation mode.
- `run_proposed_framework_comparison`: runs repeated experiments across multiple evaluation methods and collects performance statistics.

Included accuracy-analysis functions:
- `evaluate_accuracy_from_history`: performs offline comparison between local and full evaluation objectives.
- `build_accuracy_log_from_history`: constructs detailed iteration-level accuracy logs.

Included export functions:
- `export_comparison_to_excel`: exports summary tables, run statistics, accuracy results, lineplans, and iteration logs to Excel.

Compared evaluation modes include:
- `:full`
- `:local_dijkstra_pq`
- `:local_floyd_subset`
- `:threshold_basic`
- `:adaptive_local`


## Analysis tool

### `lineplan_check.jl`:
Contains a standalone lineplan feasibility and objective analysis tool. Integrates structural feasibility checking, objective evaluation, solution diagnostics, and Excel export procedures.

Main helper functions:
- `compute_present_absent_stops`: computes present and absent stop indicators for a lineplan.

Included feasibility functions:
- `explain_structure_only`: checks structural feasibility of a lineplan (route length constraints, network connectivity, and valid loop/repeated-node structure) and reports the first detected violation.

Included evaluation functions:
- `compute_objective_structure_only`: computes objective values and shortest-path information for structurally feasible lineplans.

Included diagnostic functions:
- `solution_diagnostics`: computes diagnostics related to OD connectivity, served demand, and stop coverage.

Included export procedures:
- console output of feasibility, objective values, and diagnostics.
- Excel export of evaluation and diagnostic results.
