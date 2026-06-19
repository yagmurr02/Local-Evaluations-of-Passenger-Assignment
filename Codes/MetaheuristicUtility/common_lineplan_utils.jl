# SHARED HELPER FUNCTIONS FOR METAHEURISTIC FRAMEWORKS

# Network edge check
if !isdefined(@__MODULE__, :has_edge)
    function has_edge(network, u::Int, v::Int)
        return u != v && isfinite(network[u, v])
    end
end

# Extract the nonzero stop sequence of a line
if !isdefined(@__MODULE__, :route_from_lineplan)
    function route_from_lineplan(lineplan, L::Int)
        return Int.(filter(!=(0), vec(lineplan[L, :])))
    end
end

# Extract the stop sequence of a line
if !isdefined(@__MODULE__, :line_stops)
    function line_stops(lineplan::AbstractMatrix{Int}, line_id::Int)
        return route_from_lineplan(lineplan, line_id)
    end
end

# Overwrite a line route in the line plan
if !isdefined(@__MODULE__, :set_route!)
    function set_route!(lineplan::AbstractMatrix{Int}, line_id::Int, route::Vector{Int})
        lineplan[line_id, :] .= 0
        lineplan[line_id, 1:length(route)] = route
        return lineplan
    end
end

# Check whether a route forms a proper loop line
if !isdefined(@__MODULE__, :is_proper_loop)
    function is_proper_loop(route::Vector{Int})
        n = length(route)

        if n < 3
            return false
        end

        if route[1] != route[end]
            return false
        end

        core = route[1:end-1]
        return length(unique(core)) == length(core)
    end
end

# Check whether a route satisfies node-count, uniqueness/loop, and connectivity conditions
if !isdefined(@__MODULE__, :is_valid_route)
    function is_valid_route(route::Vector{Int},
                            network;
                            min_nodes::Int = 1,
                            max_nodes::Int = typemax(Int),
                            allow_loop_lines::Bool = false)

        n = length(route)

        if n < min_nodes || n > max_nodes
            return false
        end

        if length(unique(route)) == n
            # simple path: valid
        elseif allow_loop_lines && is_proper_loop(route)
            # proper loop: valid
        else
            return false
        end

        if n >= 2
            @inbounds for i in 1:(n - 1)
                if !has_edge(network, route[i], route[i + 1])
                    return false
                end
            end
        end

        return true
    end
end