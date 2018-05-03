_NI() = error("Not implemented")

# Interface for ordering
abstract type Ordering
end

get_var(ordering::Ordering, layer::Int) = _NI


# Simple ordering implementations

# No ordering: variable corresponds to layer
struct NoOrdering <: Ordering
end

function get_var(ordering::NoOrdering, layer::Int)
    return layer
end


# Fixed ordering
struct FixedOrdering <: Ordering
    fixed_ordering::Array{Int}
end

function FixedOrdering(fixed_ordering::Array{Int})
    # Check if the ordering is a valid permutation
    if !isperm(fixed_ordering)
            error("Invalid fixed ordering")
    end
    return new(fixed_ordering)
end

"""
#Alternative definition for FixedOrdering
function FixedOrdering(fixed_ordering::Array{Int})
    # Check if the ordering is a valid permutation
    n = length(fixed_ordering)
    ordering_check = zeros(n)
    for v in fixed_ordering
        if v <= 0 || v > n
            error("Invalid fixed ordering: values outside 1:n range")
        end
        if ordering_check[v] != 0
            error("Invalid fixed ordering: not a permutation of 1:n")
        end
        ordering_check[v] = 1
    end
    return new(fixed_ordering)
end
"""

function get_var(ordering::FixedOrdering, layer::Int)
    return ordering.fixed_ordering[layer]
end
