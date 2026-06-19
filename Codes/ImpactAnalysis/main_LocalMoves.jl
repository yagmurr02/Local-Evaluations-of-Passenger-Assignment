include(joinpath(@__DIR__, "..", "CoreImplementation", "evaluator.jl"))
include(joinpath(@__DIR__, "..", "CoreImplementation", "LocalMoves_BEM.jl"))
include("RandomExperiments.jl")

using DataFrames, DataStructures, DelimitedFiles, Random, Statistics, Tables, XLSX

# Load base data
tp = 5
basepath = @__DIR__

coords  = readdlm(joinpath(basepath, "MandlCoords.txt"))
OD      = readdlm(joinpath(basepath, "MandlDemand.txt"))
network = readdlm(joinpath(basepath, "MandlNetwork.txt"))

lineplan = [
    1 2 3 6 8 10 11 12 0 0;
    1 2 5 4 6 8 10 11 0 0;
    10 7 15 6 3 2 4 12 0 0;
    9 15 8 10 14 13 11 12 0 0
]

# Compute initial DLN
shortest_old, lineinfo_old, directlink_old, path_old = shortest_paths(network, lineplan, tp)

N = 1000
Random.seed!(42)

# Insertion at Terminal
results_insertion_at_terminal = run_random_experiments_insertion_at_terminal("Insertion at Terminal", N, lineplan, network, tp, OD)
df_insertion_at_terminal = results_to_dataframe(results_insertion_at_terminal)
summary_insertion_at_terminal = summarize_results(df_insertion_at_terminal)
write_sheet!("LocalMoves_Results.xlsx", "insertion_terminal_results", df_insertion_at_terminal; overwrite_file=true)
write_sheet!("LocalMoves_Results.xlsx", "insertion_terminal_summary", summary_insertion_at_terminal; overwrite_file=false)

# Insertion
results_insertion = run_random_experiments_insertion("Insertion", N, lineplan, network, tp, OD)
df_insertion = results_to_dataframe(results_insertion)
summary_insertion = summarize_results(df_insertion)
write_sheet!("LocalMoves_Results.xlsx", "insertion_results", df_insertion; overwrite_file=false)
write_sheet!("LocalMoves_Results.xlsx", "insertion_summary", summary_insertion; overwrite_file=false)

# Removal at Terminal
results_removal_at_terminal = run_random_experiments_removal_at_terminal("Removal at Terminal", N, lineplan, network, tp, OD)
df_removal_at_terminal = results_to_dataframe(results_removal_at_terminal)
summary_removal_at_terminal = summarize_results(df_removal_at_terminal)
write_sheet!("LocalMoves_Results.xlsx", "removal_terminal_results", df_removal_at_terminal; overwrite_file=false)
write_sheet!("LocalMoves_Results.xlsx", "removal_terminal_summary", summary_removal_at_terminal; overwrite_file=false)

# Removal
results_removal = run_random_experiments_removal("Removal", N, lineplan, network, tp, OD)
df_removal = results_to_dataframe(results_removal)
summary_removal = summarize_results(df_removal)
write_sheet!("LocalMoves_Results.xlsx", "removal_results", df_removal; overwrite_file=false)
write_sheet!("LocalMoves_Results.xlsx", "removal_summary", summary_removal; overwrite_file=false)

# Swap within Line
results_swap_within_line = run_random_experiments_swap_within_line("Swap Within Line", N, lineplan, network, tp, OD)
df_swap_within_line = results_to_dataframe(results_swap_within_line)
summary_swap_within_line = summarize_results(df_swap_within_line)
write_sheet!("LocalMoves_Results.xlsx", "swap_within_results", df_swap_within_line; overwrite_file=false)
write_sheet!("LocalMoves_Results.xlsx", "swap_within_summary", summary_swap_within_line; overwrite_file=false)

# Swap Between Lines
results_swap_between_lines = run_random_experiments_swap_between_lines("Swap Between Lines", N, lineplan, network, tp, OD)
df_swap_between_lines = results_to_dataframe(results_swap_between_lines)
summary_swap_between_lines = summarize_results(df_swap_between_lines)
write_sheet!("LocalMoves_Results.xlsx", "swap_between_results", df_swap_between_lines; overwrite_file=false)
write_sheet!("LocalMoves_Results.xlsx", "swap_between_summary", summary_swap_between_lines; overwrite_file=false)

# Transfer within Line
results_transfer_within_line = run_random_experiments_transfer_within_line("Transfer Within Line", N, lineplan, network, tp, OD)
df_transfer_within_line = results_to_dataframe(results_transfer_within_line)
summary_transfer_within_line = summarize_results(df_transfer_within_line)
write_sheet!("LocalMoves_Results.xlsx", "transfer_within_results", df_transfer_within_line; overwrite_file=false)
write_sheet!("LocalMoves_Results.xlsx", "transfer_within_summary", summary_transfer_within_line; overwrite_file=false)

# Transfer between Lines
results_transfer_between_lines = run_random_experiments_transfer_between_lines("Transfer Between Lines", N, lineplan, network, tp, OD)
df_transfer_between_lines = results_to_dataframe(results_transfer_between_lines)
summary_transfer_between_lines = summarize_results(df_transfer_between_lines)
write_sheet!("LocalMoves_Results.xlsx", "transfer_between_results", df_transfer_between_lines; overwrite_file=false)
write_sheet!("LocalMoves_Results.xlsx", "transfer_between_summary", summary_transfer_between_lines; overwrite_file=false)

# Substitution
results_substitution = run_random_experiments_substitution("Substitution", N, lineplan, network, tp, OD)
df_substitution = results_to_dataframe(results_substitution)
summary_substitution = summarize_results(df_substitution)
write_sheet!("LocalMoves_Results.xlsx", "substitution_results", df_substitution; overwrite_file=false)
write_sheet!("LocalMoves_Results.xlsx", "substitution_summary", summary_substitution; overwrite_file=false)

# Partial Segment Swap
results_partial_segment_swap = run_random_experiments_partial_segment_swap("Partial Segment Swap", N, lineplan, network, tp, OD)
df_partial_segment_swap = results_to_dataframe(results_partial_segment_swap)
summary_partial_segment_swap = summarize_results(df_partial_segment_swap)
write_sheet!("LocalMoves_Results.xlsx", "partial_segment_swap_results", df_partial_segment_swap; overwrite_file=false)
write_sheet!("LocalMoves_Results.xlsx", "partial_segment_swap_summary", summary_partial_segment_swap; overwrite_file=false)

# Reversal
results_reversal = run_random_experiments_reversal("Reversal", N, lineplan, network, tp, OD)
df_reversal = results_to_dataframe(results_reversal)
summary_reversal = summarize_results(df_reversal)
write_sheet!("LocalMoves_Results.xlsx", "reversal_results", df_reversal; overwrite_file=false)
write_sheet!("LocalMoves_Results.xlsx", "reversal_summary", summary_reversal; overwrite_file=false)

# Segment Reversal
results_segment_reversal = run_random_experiments_segment_reversal("Segment Reversal", N, lineplan, network, tp, OD)
df_segment_reversal = results_to_dataframe(results_segment_reversal)
summary_segment_reversal = summarize_results(df_segment_reversal)
write_sheet!("LocalMoves_Results.xlsx", "segment_reversal_results", df_segment_reversal; overwrite_file=false)
write_sheet!("LocalMoves_Results.xlsx", "segment_reversal_summary", summary_segment_reversal; overwrite_file=false)

# Line Addition
results_line_addition = run_random_experiments_line_addition("Line Addition", N, lineplan, network, tp, OD)
df_line_addition = results_to_dataframe(results_line_addition)
summary_line_addition = summarize_results(df_line_addition)
write_sheet!("LocalMoves_Results.xlsx", "line_addition_results", df_line_addition; overwrite_file=false)
write_sheet!("LocalMoves_Results.xlsx", "line_addition_summary", summary_line_addition; overwrite_file=false)

# Line Removal
results_line_removal = run_random_experiments_line_removal("Line Removal", N, lineplan, network, tp, OD)
df_line_removal = results_to_dataframe(results_line_removal)
summary_line_removal = summarize_results(df_line_removal)
write_sheet!("LocalMoves_Results.xlsx", "line_removal_results", df_line_removal; overwrite_file = false)
write_sheet!("LocalMoves_Results.xlsx", "line_removal_summary", summary_line_removal; overwrite_file = false)