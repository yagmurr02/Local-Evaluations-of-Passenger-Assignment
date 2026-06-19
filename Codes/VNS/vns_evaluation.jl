"""
VNS (Variable Neighborhood Search) evaluation module.
Uses the unified evaluation framework for both SA and VNS moves.
"""

include(joinpath(@__DIR__, "vns_moves.jl"))
include(joinpath(@__DIR__, "evaluation.jl"))
