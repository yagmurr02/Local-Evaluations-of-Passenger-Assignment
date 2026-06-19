if !isdefined(@__MODULE__, :ProposedMove)

abstract type ProposedMove end

struct ProposedInsertion <: ProposedMove
    line::Int
    left_node::Int
    right_node::Int
    inserted_node::Int
end

struct ProposedRemoval <: ProposedMove
    line::Int
    removed_node::Int
end

struct ProposedSwapWithinLine <: ProposedMove
    line::Int
    node1::Int
    node2::Int
end

struct ProposedSubstitution <: ProposedMove
    line::Int
    old_node::Int
    new_node::Int
end

struct ProposedReversal <: ProposedMove
    line::Int
end

struct ProposedInsertionTerminal <: ProposedMove
    line::Int
    node::Int
    position::String
end

struct ProposedRemovalTerminal <: ProposedMove
    line::Int
    node::Int
    position::String
end

struct ProposedSegmentReversal <: ProposedMove
    line::Int
    start_node::Int
    end_node::Int
end

struct ProposedSwapBetweenLines <: ProposedMove
    line1::Int
    node1::Int
    line2::Int
    node2::Int
end

struct ProposedPositionChangeBetweenLines <: ProposedMove
    from_line::Int
    to_line::Int
    moved_node::Int
    left_node::Int
    right_node::Int
end

struct ProposedPartialLineSwap <: ProposedMove
    line1::Int
    node1a::Int
    node1b::Int
    line2::Int
    node2a::Int
    node2b::Int
end

end