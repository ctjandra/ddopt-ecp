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

function get_var(ordering::FixedOrdering, layer::Int)
    return ordering.fixed_ordering[layer]
end
