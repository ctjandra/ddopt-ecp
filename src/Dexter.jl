__precompile__()

module Dexter
    import MathProgBase

    include("solver.jl")

    include("dd/ordering.jl")
    include("dd/construction.jl")
    include("dd/modification.jl")
    include("dd/computation.jl")

    include("dd/mathprog/decomposition.jl")
    include("dd/mathprog/dp_specification.jl")

    include("ip/oa/ecp.jl")
end
