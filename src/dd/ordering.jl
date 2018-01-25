_NI() = error("Not implemented")

# Interface for ordering
abstract type Ordering
end

get_var(ordering::Ordering, layer::Int) = _NI


# Simple ordering implementations

# TODO Implement orderings

# No ordering: variable corresponds to layer
struct NoOrdering <: Ordering
end

function get_var(ordering::NoOrdering, layer::Int)
    return layer
end
