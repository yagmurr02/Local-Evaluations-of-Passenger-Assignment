abstract type VNSMove end

struct VNSInsertion <: VNSMove
    line::Int
    left_node::Int
    right_node::Int
    inserted_node::Int
end

struct VNSRemoval <: VNSMove
    line::Int
    removed_node::Int
end

struct VNSSwapWithinLine <: VNSMove
    line::Int
    node1::Int
    node2::Int
end

struct VNSSubstitution <: VNSMove
    line::Int
    old_node::Int
    new_node::Int
end

struct VNSPartialSegmentSwap <: VNSMove
    line1::Int
    node1a::Int
    node1b::Int
    line2::Int
    node2a::Int
    node2b::Int
end

struct VNSReversal <: VNSMove
    line::Int
end

struct VNSInsertionTerminal <: VNSMove
    line::Int
    node::Int
    position::String   # "start" or "end"
end

struct VNSRemovalTerminal <: VNSMove
    line::Int
    node::Int
    position::String   # "start" or "end"
end