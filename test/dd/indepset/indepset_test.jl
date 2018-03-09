include("../../../src/Dexter.jl")
using Dexter

using Gadfly
using GraphPlot

include("read_graph.jl")
include("indepset.jl")

# Generate graph
# g = read_graph("test/dd/indepset/instances/random_200_50_0.clq")
# g = PathGraph(10)
# g = erdos_renyi(20, 0.5)
g = watts_strogatz(20, 7, 0.8)

# Plot graph
gplot(g) |> PNG("graph.png", 16cm, 16cm)

# Create independent set specs for decision diagram
initial_state = create_indepset_initial_state(g)
transition = create_indepset_transition_function(g)
domain_function = Dexter.domain_function_binary
indepset_specs = Dexter.ProblemSpecs(initial_state, transition, domain_function)

# Set ordering and objective
n = nv(g)
ordering = Dexter.NoOrdering()
objective = ones(n)

# Construct decision diagram
dd = Dexter.construct_DD(n, indepset_specs, ordering=ordering, objective=objective)

# Plot decision diagram
locs_x = [Dexter.get_node_idlayer(dd, i) * 1.0 for i in 1:nv(dd.graph)]
locs_y = [Dexter.get_node_layer(dd, i) / n for i in 1:nv(dd.graph)]
gplot(dd.graph, locs_x, locs_y) |> PNG("test.png", 16cm, 16cm)
