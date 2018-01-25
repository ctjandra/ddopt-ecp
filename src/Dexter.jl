__precompile__()

module Dexter
    import MathProgBase

    include("solver.jl")

    include("dd/ordering.jl")
    include("dd/construction.jl")
end
