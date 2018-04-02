__precompile__()

module Dexter
    import MathProgBase

    include("solver.jl")

    include("dd/ordering.jl")
    include("dd/construction.jl")
    include("dd/functions.jl")

    include("cut/cutgenerator.jl")
end
