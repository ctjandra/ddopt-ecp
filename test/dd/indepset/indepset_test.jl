using Gadfly
using GraphPlot

include("read_graph.jl")

g = read_graph("test/dd/indepset/instances/random_200_50_0.clq")

gplot(g) |> PNG("test.png", 16cm, 16cm)
