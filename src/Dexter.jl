__precompile__()

module Dexter
    import MathProgBase

    include("solver.jl")

    include("dd/ordering.jl")
    include("dd/core.jl")
    include("dd/construction.jl")
    include("dd/computation.jl")

    include("dd/bb/bb_tree.jl")
    include("dd/mathprog/decomposition.jl")
    include("dd/mathprog/function_analysis.jl")
    include("ip/oa/ecp.jl")
end
