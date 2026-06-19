# Extract nonzero stops from a line as a clean vector
get_stops(lineplan, line_id) = filter(!=(0), lineplan[line_id, :])

# Write a vector of stops back into the lineplan (zero-padded)
function set_stops!(lineplan, line_id, stops)
    MAX = size(lineplan, 2)

    # Truncate if segment is too long
    if length(stops) > MAX
        stops = stops[1:MAX]
    end

    # Write real nodes
    lineplan[line_id, 1:length(stops)] .= stops

    # Fill the rest with zeros
    if length(stops) < MAX
        lineplan[line_id, length(stops)+1:MAX] .= 0
    end
end

# Insert(line, node, position): 
function Insert_node(stops::Vector{Int}, node::Int, position::Int)
    # insert before position
    return vcat(stops[1:position-1], node, stops[position:end])
end

# Remove(line, node): 
function Remove_node(stops::Vector{Int}, node::Int)
    return filter(!=(node), stops)
end

# InsertSegment(line, segment, after node): 
function InsertSegment_after(stops::Vector{Int}, segment::Vector{Int}, node::Int)
    idx = findfirst(==(node), stops)
    idx === nothing && error("Node $node not found for InsertSegment.")

    return vcat(
        stops[1:idx],          # everything up to the anchor
        segment,               # inserted segment
        stops[idx+1:end]       # remainder
    )
end

# RemoveSegment(line, start node, end node): 
function RemoveSegment_between(stops::Vector{Int}, node_i::Int, node_j::Int)
    idx_i = findfirst(==(node_i), stops)
    idx_j = findfirst(==(node_j), stops)

    (idx_i === nothing || idx_j === nothing) &&
        error("Nodes $node_i or $node_j not found in stops.")

    # ensure node_i comes before node_j
    if idx_i > idx_j
        idx_i, idx_j = idx_j, idx_i
    end

    # nothing between them
    if idx_i + 1 > idx_j - 1
        return stops, Int[]
    end

    segment   = stops[idx_i+1 : idx_j-1]            # middle part
    new_stops = vcat(stops[1:idx_i], stops[idx_j:end])  # keep i and j, remove middle

    return new_stops, segment
end

#  Reverse(segment): 
Reverse_segment(segment::Vector{Int}) = reverse(segment)

# CreateSegment([node,...,node]):
CreateSegment(nodes::Vector{Int}) = copy(nodes)