include(joinpath(@__DIR__, "common_lineplan_utils.jl"))

abstract type Move end

struct TerminalInsertion <: Move
    line::Int
    node::Int
end

struct TerminalRemoval <: Move
    line::Int
end

struct Reversal <: Move
    line::Int
end

# Apply terminal insertion
function apply_move(lineplan, move::TerminalInsertion)
    new_lineplan = deepcopy(lineplan)
    route = route_from_lineplan(new_lineplan, move.line)

    push!(route, move.node)

    set_route!(new_lineplan, move.line, route)
    return new_lineplan
end

# Apply terminal removal
function apply_move(lineplan, move::TerminalRemoval)
    new_lineplan = deepcopy(lineplan)
    route = route_from_lineplan(new_lineplan, move.line)

    if length(route) <= 1
        return nothing
    end

    route = route[2:end]

    set_route!(new_lineplan, move.line, route)
    return new_lineplan
end

# Apply reversal
function apply_move(lineplan, move::Reversal)
    new_lineplan = deepcopy(lineplan)
    route = route_from_lineplan(new_lineplan, move.line)

    reverse!(route)

    set_route!(new_lineplan, move.line, route)
    return new_lineplan
end